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

/**
 * \file
 * A header-only PID controller based largely on <a
 * href="http://www.cds.caltech.edu/~murray/amwiki/index.php/PID_Control">
 * Chapter 10</a> of <a href="http://www.worldcat.org/isbn/0691135762">Astrom
 * and Murray</a>.
 *
 * This proportional-integral-derivative (PID) controller features
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
 * \image html  helm.png "Block diagram for the controller"
 * \image latex helm.eps "Block diagram for the controller" width=\textwidth
 *
 * Let \f$f\f$ be a first-order, low-pass filtered version of controlled process
 * output \f$y\f$ governed by
 * \f{align}{
 *     \frac{\mathrm{d}}{\mathrm{d}t} f &= \frac{y - f}{T_f}
 * \f}
 * where \f$T_f\f$ is a filter time scale.  Then, in the time domain and
 * expressed in positional form, the control signal \f$v\f$ evolves according to
 * \f{align}{
 *     v(t) &= k_p \, e(t)
 *           + k_i \int_0^t e(t) \,\mathrm{d}t
 *           + k_t \int_0^t e_s(t) \,\mathrm{d}t
 *           - k_d \frac{\mathrm{d}}{\mathrm{d}t} f(t)
 *     \\
 *          &= k_p \left[
 *                 e(t)
 *               + \frac{1}{T_i} \int_0^t e(t) \,\mathrm{d}t
 *               + \frac{1}{T_t} \int_0^t e_s(t) \,\mathrm{d}t
 *               - T_d \frac{\mathrm{d}}{\mathrm{d}t} f(t)
 *             \right]
 *     \\
 *          &= k_p \left[
 *                 \left(r(t) - y(t)\right)
 *               + \frac{1}{T_i} \int_0^t \left(r(t) - y(t)\right) \,\mathrm{d}t
 *               + \frac{1}{T_t} \int_0^t \left(u(t) - v(t)\right) \,\mathrm{d}t
 *               + \frac{T_d}{T_f}\left(f(t) - y(t)\right)
 *             \right]
 * \f}
 * where \f$u\f$ is the actuator position and \f$r\f$ is the desired reference
 * or "setpoint" value.  Constants \f$T_i\f$, \f$T_t\f$, and \f$T_d\f$ are the
 * integral, automatic reset, and derivative time scales while \f$k_p\f$
 * specifies the unified gain.  Differentiating one finds the "incremental" form
 * written for continuous time,
 * \f{align}{
 *     \frac{\mathrm{d}}{\mathrm{d}t} v(t) &= k_p \left[
 *               - \frac{\mathrm{d}}{\mathrm{d}t} y(t)
 *               + \frac{r(t) - y(t)}{T_i}
 *               + \frac{u(t) - v(t)}{T_t}
 *               + \frac{T_d}{T_f}\left(
 *                   \frac{\mathrm{d}}{\mathrm{d}t} f(t)
 *                 - \frac{\mathrm{d}}{\mathrm{d}t} y(t)
 *                 \right)
 *             \right]
 *             .
 * \f}
 * Here, to avoid controller kick on instantaneous reference value changes, we
 * assume \f$\frac{\mathrm{d}}{\mathrm{d}t} r(t) = 0\f$.  This assumption is
 * sometimes called "derivative on measurement" in reference to neglecting the
 * non-measured portion of the error derivative
 * \f$\frac{\mathrm{d}}{\mathrm{d}t} e(t)\f$.
 *
 * Obtaining a discrete time evoluation equation is straightforward.  Multiply
 * the above continuous result by the time differential, substitute first
 * order backward differences, and incorporate the low-pass filter in a
 * consistent fashion.  One then finds the following:
 * \f{align}{
 *     {\mathrm{d}t}_i &= t_i - t_{i-1}
 * \f}
 * \f{align}{
 *     f(t_i) &= \frac{ {\mathrm{d}t}_i\,y(t_i) + T_f\,f(t_{i-1}) }
 *                 { T_f + {\mathrm{d}t}_i }
 *             = \alpha y(t_i) + (1 - \alpha) f_{i-1}
 *               \quad\text{with }
 *               \alpha=\frac{{\mathrm{d}t}_i}{T_f + {\mathrm{d}t}_i}
 * \f}
 * \f{align}{
 *     {\mathrm{d}f}_i &= f(t_i) - f(t_{i-1})
 *                      = \alpha\left( y(t_i) - f(t_{i-1}) \right)
 * \f}
 * \f{align}{
 *     {\mathrm{d}y}_i &= y(t_i) - y(t_{i-1})
 * \f}
 * \f{align}{
 *     {\mathrm{d}v}_i &= k_p \left[
 *                 {\mathrm{d}t}_i \left(
 *                   \frac{r(t_i) - y(t_i)}{T_i}
 *                 + \frac{u(t_i) - v(t_i)}{T_t}
 *                 \right)
 *               + \frac{T_d}{T_f}\left(
 *                   {\mathrm{d}f}_i - {\mathrm{d}y}_i
 *                 \right)
 *               - {\mathrm{d}y}_i
 *             \right]
 * \f}
 * where notice \f$f(t)\f$ is nothing but an exponential weighted moving average
 * of \f$y(t)\f$ that permits varying the sampling rate.  An implementation
 * needs only to track two pieces of state, namely \f$f(t_{i-1})\f$ and
 * \f$y(t_{i-1})\f$, across time steps.
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
 * Tuning parameters and internal state for an incremental PID controller.
 *
 * Gain #kp has units of \f$u_0 / y_0\f$ where \f$u_0\f$ and \f$y_0\f$
 * are the natural actuator and process observable signals, respectively.
 * Parameter #Tt has units of time multiplied by \f$u_0 / y_0\f$.
 * Parameters #Td, #Tf, and #Ti possess units of time.  Time units are
 * fixed by the scaling provided in the \c dt argument to helm_steady().
 */
struct helm_state
{
    double kp;  /**< Proportional gain modifying P, I, and D terms.    */
    double Td;  /**< Time scale governing derivative action.
                     Set to zero to disable derivative control.        */
    double Tf;  /**< Time scale filtering process observable for D.
                     Set to infinity to disable observable filtering.  */
    double Ti;  /**< Time scale governing integral action.
                     Set to infinity to disable integral control.      */
    double Tt;  /**< Time scale governing automatic reset.
                     Set to infinity to disable automatic reset.       */
    double y;   /**< Internal tracking of the process observable.      */
    double f;   /**< Internal tracking the filtered process.           */
};

/**
 * \brief Reset all tuning parameters, but \e not transient state.
 *
 * Resets gain to one and disables filtering, integral action, and derivative
 * action.  Enable those terms by setting their associated time scales.
 *
 * \param[in,out] h Houses tuning parameters to be reset.
 * \return Argument \c h to permit call chaining.
 */
static inline
struct helm_state *
helm_reset(struct helm_state * const h)
{
    h->kp = 1;        // Unit gain
    h->Td = 0;        // No derivative action
    h->Tf = INFINITY; // No filtering
    h->Ti = INFINITY; // No integral action
    h->Tt = INFINITY; // No automatic reset
    return h;
}

/**
 * \brief Reset any transient state, but \e not tuning parameters.
 *
 * Necessary to achieve bumpless manual-to-automatic transitions
 * before calling to helm_steady() after a period of manual control,
 * including \e before the first call to helm_steady().
 *
 * \param[in,out] h Houses transient state to be reset.
 * \return Argument \c h to permit call chaining.
 */
static inline
struct helm_state *
helm_approach(struct helm_state * const h)
{
    assert(h->Td >= 0);
    assert(h->Tf >  0);
    assert(h->Ti >  0);
    assert(h->Tt >  0);
    h->f = NAN;
    return h;
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
 * \see Overview of \ref helm.h for the discrete evolution equations.
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
    double dv = 0;

    if (!isnan(y)) {                      // Avoid driving blind

        if (isnan(h->f)) {                // Avoid startup kick
            h->y = y;
            h->f = y;
        }

        double a, df, dy;
        a   = dt / (h->Tf + dt);          // Convex combination parameter alpha
        df  = a*(y - h->f);               // Filtered difference for y
        dy  =    y - h->y ;               // Backward difference for y
        dv += (r - y) / h->Ti;            // Action from integral control
        dv += (u - v) / h->Tt;            // Action from automatic reset
        dv *= dt;                         // Scale integral actions by time step
        dv += (h->Td / h->Tf)*(df - dy);  // Action from derivative control
        dv += /*dr=0*/ - dy;              // Action from proporational control
        dv *= h->kp;                      // Scale by unified gain parameter

        h->y  = y;                        // Update observable for next call
        h->f += df;                       // Update filter for next call
    }

    return dv;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HELM_H */
