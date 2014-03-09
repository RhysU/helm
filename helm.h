//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef HELM_H
#define HELM_H

#include <assert.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

// TODO Document the continuous time and discretized equations
// FIXME Check claims on no filtering with Tf = 1

#ifndef HELM_SUPPRESS_DOXYGEN_MAINPAGE
/** \mainpage
 * Please see \ref helm.h for design details and http://github.com/RhysU/helm
 * for project information.
 */
#endif

/**
 * \file
 * \brief A header-only C99 proportional-integral-derivative (PID) controller.
 *
 * \image html  helm.png "Controller block diagram"
 * \image latex helm.eps "Controller block diagram" width=\textwidth
 * The controller features
 * <ul>
 *   <li>low pass filtering of the process derivative,</li>
 *   <li>windup protection,</li>
 *   <li>automatic reset on actuator saturation,</li>
 *   <li>anti-kick on setpoint change using "derivative on measurement",</li>
 *   <li>incremental output for bumpless manual-to-automatic transitions,</li>
 *   <li>a unified controller gain parameter,</li>
 *   <li>exposure of all independent physical time scales, and</li>
 *   <li>the ability to accommodate varying sample rate.</li>
 * </ul>
 * The design and nomenclature is based largely on Figure 10.11 in <a
 * href="http://www.worldcat.org/isbn/9781400828739">Feedback Systems</a> by
 * Astrom and Murray.
 *
 * In the time domain, the governing positional equation is
 * \f{align}{
 *  v(t) &= k_p e(t)
 *        + k_i \int_0^t e(t) \mathrm{d}t
 *        + k_t \int_0^t e_s(t) \mathrm{d}t
 *        - k_d \frac{\mathrm{d}}{\mathrm{d}t} f(t)
 * \f}
 * where
 * \f{align}{
 *  T_f \frac{\mathrm{d}}{\mathrm{d}t} f &= y - f
 * \f}
 * provides low order filtering.
 *
 * Sample written with nomenclature from helm_state() and helm_steady():
 * \code
 *   struct helm_state h;
 *
 *   // Set PID parameters from commonly given \c kp, \c ki, \c kt, and \c kd
 *   helm_reset(h);
 *   h->kp = kp;
 *   h->Td = kd / h->kp;
 *   h->Tf = h->Td / 10;  // Astrom and Murray p.308 suggests 2--20
 *   h->Ti = h->kp / ki;
 *   h->Tt = h->kp / kt;
 *
 *   // Enable automatic control and evolve
 *   helm_approach(h);
 *   for (int i = 0; i < N; ++i) {
 *      y  = process(dt, u);
 *      v += helm_steady(h, dt, r, u, v, y);
 *      u  = actuate(dt, v);
 *   }
 *
 *   // Disable controller and evolve
 *   for (int i = 0; i < N; ++i) {  // E
 *      y  = process(dt, u);
 *      u  = actuate(dt, v);
 *   }
 *
 *   // Re-enable automatic control and evolve
 *   helm_approach(h);
 *   for (int i = 0; i < N; ++i) {
 *      y  = process(dt, u);
 *      v += helm_steady(h, dt, r, u, v, y);
 *      u  = actuate(dt, v);
 *   }
 * \endcode
 */

/**
 * State for an incremental PID controller, including all tuning parameters.
 */
struct helm_state
{
    /**
     * Controller tuning parameters.
     *
     * Gain has units of <code>u0 / y0</code>.
     * Tt has units of time multiplied by <code>u0 / y0</code>.
     * All other time scales possess units of time.
     *
     * Setting a time scale to \c INFINITY disables the associated term.
     * 
     * @{
     */
    double kp;  /**< Proportional gain modifying P, I, and D terms.  */
    double Td;  /**< Time scale governing derivative action.         */
    double Tf;  /**< Time scale filtering process observable for D.  */
    double Ti;  /**< Time scale governing integral action.           */
    double Tt;  /**< Time scale governing automatic reset.           */
    /**@}*/

    /**
     * Internal state maintained between calls to helm_steady().
     * @{
     */
    double y;   /**< Tracks instantaneous process observable. */
    double f;   /**< Tracks filtered process observable.      */
    /**@}*/
};

/**
 * \brief Reset all tuning parameters, but \e not transient state.
 *
 * Resets gain to one and  disables filtering, integral action, and derivative
 * action.  Enable those terms by setting the associated time scales.
 */
static inline
void
helm_reset(struct helm_state * const h)
{
    h->kp = 1;        // Unit gain
    h->Td = INFINITY; // No derivative action
    h->Tf = 1;        // No filtering
    h->Ti = INFINITY; // No integral action
    h->Tt = INFINITY; // No automatic reset
}

/**
 * \brief Forget any transient state, but \e not tuning parameters.
 *
 * Necessary to achieve bumpless manual-to-automatic transitions
 * before calling to helm_steady() after a period of manual control.
 */
static inline
void
helm_approach(struct helm_state * const h)
{
    assert(h->Td > 0);
    assert(h->Tf > 0);
    assert(h->Ti > 0);
    assert(h->Tt > 0);
    h->f = NAN;
}

/**
 * \brief Find the control signal necessary to steady unsteady process y(t).
 *
 * \param[in,out] h  Tuning parameters and state maintained across invocations.
 * \param[in]     dt Time since last samples collected.
 * \param[in]     r  Reference value, often called the "setpoint".
 * \param[in]     u  Actuator signal currently observed.
 * \param[in]     v  Actuator signal currently requested.
 * \param[in]     y  Observed process output to drive to \c r.
 *
 * \return Incremental suggested change to control signal \c v.
 */
static inline
double
helm_steady(struct helm_state * const h,
            const double dt,
            const double r,
            const double u,
            const double v,
            const double y)
{
    double dy, df, dv = 0;

    if (!isnan(y)) {                      // Avoid driving blind

        if (isnan(h->f)) {                // Avoid startup kick
            h->y = y;
            h->f = y;
        }

        dy  = y - h->y;                   // Backward difference for y
        df  = (dt / h->Tf)*(y - h->f);    // Filtered difference for y
        dv += (r - y) / h->Ti;            // Action from integral control
        dv += (u - v) / h->Tt;            // Action from automatic reset
        dv *= dt;                         // Scale integral actions by time step
        dv += (h->Td / h->Tf)*(df - dy);  // Action from derivative control
        dv += /*dr=0*/ - dy;              // Action from proporational control
        dv *= h->kp;                      // Scale by unified gain parameter

        h->y  = y;                        // Track observable for next call
        h->f += df;                       // Track filter for next call
    }

    return dv;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HELM_H */
