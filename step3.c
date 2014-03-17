//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** \file
 * Sample controlling a 3rd order system across a unit step in the setpoint.
 *
 * This sample can be used to test controller behavior against known good
 * results.  For example, those presented in Figure 10.2 within <a
 * href="http://www.cds.caltech.edu/~murray/amwiki/index.php/PID_Control">
 * Chapter 10</a> of <a href="http://www.worldcat.org/isbn/0691135762">Astrom
 * and Murray</a>.
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "helm.h"

/**
 * Advance the temporal state of a model given by transfer function \f$
 * \frac{y(s)}{u(s)} = \frac{b_0}{s^3 + a_2 s^2 + a_1 s + a_0} \f$.
 *
 * Given the process transfer function
 * \f{align}{
 *     \frac{y(s)}{u(s)} = \frac{b_0}{s^3 + a_2 s^2 + a_1 s + a_0}
 * \f}
 * a matching state space model consisting of 1st order differential equations
 * <a
 * href="http://lpsa.swarthmore.edu/Representations/SysRepTransformations/SysRepTransfAll.html">
 * can be derived</a> with the form
 * \f{align}{
 *     \frac{\mathrm{d}}{\mathrm{d}t}
 *     \begin{bmatrix} y_0(t) \\ y_1(t) \\ y_2(t) \end{bmatrix}
 *     &=
 *     \begin{bmatrix}
 *            0 &    1 &    0 \\
 *            0 &    0 &    1 \\
 *         -a_0 & -a_1 & -a_2
 *     \end{bmatrix}
 *     \begin{bmatrix} y_0(t) \\ y_1(t) \\ y_2(t) \end{bmatrix}
 *     +
 *     \begin{bmatrix}   0 \\   0 \\ b_0 \end{bmatrix} u(t)
 * \f}
 * for constants \f$a_0\f$, \f$a_1\f$, \f$a_2\f$, and \f$b_0\f$ and time-varying
 * input data \f$u(t)\f$.  Using a semi-implicit Euler integration scheme,
 * \f{align}{
 *     \vec{y}\left(t_{i+1}\right)
 *     &=
 *       \vec{y}\left(t_i\right)
 *     + h \vec{f}\left(\vec{y}\left(t_{i+1}\right), u\left(t_i\right)\right)
 *     ,
 * \f}
 * yields a constant-coefficient linear problem for advancing by time \f$h\f$:
 * \f{align}{
 *     \begin{bmatrix}
 *            1 &   -h &      0 \\
 *            0 &    1 &     -h \\
 *         ha_0 & ha_1 & 1+ha_2
 *     \end{bmatrix}
 *     \begin{bmatrix} y_0(t_{i+1})\\y_1(t_{i+1})\\y_2(t_{i+1}) \end{bmatrix}
 *     &=
 *     \begin{bmatrix} y_0(t_i) \\ y_1(t_i) \\ y_2(t_i) \end{bmatrix}
 *     +
 *     \begin{bmatrix}   0 \\   0 \\ b_0 \end{bmatrix} u(t_i)
 * \f}
 * Left multiplying by the matrix cofactor and dividing by the determinant gives
 * a form amenable to computation,
 * \f{align}{
 *     \begin{bmatrix} y_0(t_{i+1})\\y_1(t_{i+1})\\y_2(t_{i+1}) \end{bmatrix}
 *     &=
 *     \frac{
 *       \begin{bmatrix}
 *        h (a_2+a_1 h)+1 & h (a_2 h+1)    & h^2 \\
 *        -a_0 h^2        & a_2 h+1        & h \\
 *        -a_0 h          & -h (a_1+a_0 h) & 1
 *       \end{bmatrix}
 *       \left(
 *           \begin{bmatrix} y_0(t_i) \\ y_1(t_i) \\ y_2(t_i) \end{bmatrix}
 *           +
 *           \begin{bmatrix}   0 \\   0 \\ b_0 \end{bmatrix} u(t_i)
 *       \right)
 *     }{h (h (a_0 h+a_1)+a_2)+1}
 *     .
 * \f}
 * This routine advances \f$\vec{y}(t)\f$ to \f$\vec{y}(t + h)\f$ using
 * the above result.
 *
 * @param[in    ] h Time step \f$h\f$ to be taken.
 * @param[in    ] a Coefficients \f$a_0\f$, \f$a_1\f$, and \f$a_2\f$.
 * @param[in    ] b Coefficient \f$b_0\f$.
 * @param[in    ] u Input \f$u(t)\f$.
 * @param[in,out] y On input,  state \f$y_0(t)\f$, \f$y_1(t)\f$,
 *                  and \f$y_2(t)\f$.  On output, state \f$y_0(t+h)\f$,
 *                  \f$y_1(t+h)\f$, and \f$y_2(t+h)\f$.
 */
static
void
advance(const double h,
        const double a[restrict static 3],
        const double b[restrict static 1],
        const double u[restrict static 1],
              double y[restrict static 3])
{
    const double rhs[3] = {
        y[0],
        y[1],
        y[2] + b[0]*u[0]
    };
    const double cof[3][3] = {
        { h*(a[2] + a[1]*h) + 1 , h*(a[2]*h + 1)     , h*h },
        { -a[0]*h*h             , a[2]*h + 1         , h   },
        { -a[0]*h               , -h*(a[1] + a[0]*h) , 1   }
    };
    const double det = h*(h*(a[0]*h + a[1]) + a[2]) + 1;

    for (int i = 0; i < 3; ++i) {
        y[i] = 0;
        for (int j = 0; j < 3; ++j) {
            y[i] += cof[i][j]*rhs[j];
        }
        y[i] /= det;
    }
}

static const double default_a[3] = {1, 3, 3}; ///< Default process parameters
static const double default_b[1] = {1};       ///< Default process parameters
static const double default_f    = 0.01;      ///< Default filter time scale
static const double default_kd   = 1;         ///< Default derivative gain
static const double default_ki   = 1;         ///< Default integration gain
static const double default_kp   = 1;         ///< Default proportional gain
static const double default_r    = 1;         ///< Default reference value
static const double default_t    = 1;         ///< Default time step size
static const double default_T    = 25;        ///< Default final time

/** Print usage on the given stream. */
static
void
print_usage(const char *arg0, FILE *out)
{
    fprintf(out, "Usage: %s [OPTION...]\n", arg0);
    fprintf(out, "Control 3rd-order system across a setpoint step change.\n");
    fprintf(out, "Output is tab-delimited t, u, y[0], y[1], y[2].\n");
    fputc('\n', out);
    fprintf(out, "Process transfer function"
                    "y(s)/u(s) = b0 / (s^3 + a2 s^2 + a1 s + a0):\n");
    fprintf(out, "  -0 a0\t\tSet coefficient a0 (default %g)\n", default_a[0]);
    fprintf(out, "  -1 a1\t\tSet coefficient a1 (default %g)\n", default_a[1]);
    fprintf(out, "  -2 a2\t\tSet coefficient a2 (default %g)\n", default_a[2]);
    fprintf(out, "  -b b0\t\tSet coefficient b0 (default %g)\n", default_b[0]);
    fputc('\n', out);
    fprintf(out, "Term-by-term, parallel-form PID settings\n");
    fprintf(out, "  -p kp\t\tProportional gain  (default %g)\n", default_kp);
    fprintf(out, "  -i ki\t\tIntegral gain      (default %g)\n", default_ki);
    fprintf(out, "  -d kd\t\tDerivative gain    (default %g)\n", default_kd);
    fprintf(out, "  -f Tf\t\tFilter time scale  (default %g)\n", default_f);
    fprintf(out, "  -r ref\t\tReference value   (default %g)\n", default_r);
    fputc('\n', out);
    fprintf(out, "Miscellaneous:\n");
    fprintf(out, "  -r r \t\tAdjust setpoint    (default %g)\n", default_r);
    fprintf(out, "  -t dt\t\tSet time step size (default %g)\n", default_t);
    fprintf(out, "  -T Tf\t\tSet final time     (default %g)\n", default_T);
    fprintf(out, "  -h\t\tDisplay this help and exit\n");
}

/**
 * Control the process with transfer function \f$ \frac{y(s)}{u(s)} =
 * \frac{b_0}{s^3 + a_2 s^2 + a_1 s + a_0} \f$ across a unit step change in
 * setpoint value.  That is, just prior to time zero process state \f$y(t)\f$,
 * reference value \f$r(t)\f$, actuator signal \f$u(t)\f$, and all of their
 * derivatives are zero.  At time zero, step change \f$r(t) = 1\f$ is
 * introduced.  The transfer function, in conjunction with the controller
 * dynamics, determines the controlled system response.
 */
int
main (int argc, char *argv[])
{
    // Establish mutable settings
    double a[3] = {default_a[0], default_a[1], default_a[2]};
    double b[1] = {default_b[0]};
    double f    = default_f;
    double kd   = default_kd;
    double ki   = default_ki;
    double kp   = default_kp;
    double r    = default_r;
    double t    = default_t;
    double T    = default_T;

    // Process incoming arguments
    static const char optstring[] = "0:1:2:b:d:i:p:r:t:T:h";
    for (int opt; -1 != (opt = getopt(argc, argv, optstring));) {
        switch (opt) {
        case '0': a[0] = atof(optarg);          break;
        case '1': a[1] = atof(optarg);          break;
        case '2': a[2] = atof(optarg);          break;
        case 'b': b[0] = atof(optarg);          break;
        case 'd': kd   = atof(optarg);          break;
        case 'i': ki   = atof(optarg);          break;
        case 'p': kp   = atof(optarg);          break;
        case 'r': r    = atof(optarg);          break;
        case 't': t    = atof(optarg);          break;
        case 'T': T    = atof(optarg);          break;
        case 'h': print_usage(argv[0], stdout); return EXIT_SUCCESS;
        default:  print_usage(argv[0], stderr); return EXIT_FAILURE;
        }
    }

    // Avoid infinite loops by sanitizing inputs
    if (t <= 0) {
        fprintf(stderr, "Step size t must be strictly positive\n");
        return EXIT_FAILURE;
    }
    if (T <= 0) {
        fprintf(stderr, "Final time T must be strictly positive\n");
        return EXIT_FAILURE;
    }

    // Initialize state
    double u[1] = {r};        // Actuator signal
    double v[1] = {r};        // Control signal
    double y[3] = {0, 0, 0};  // Model state

    // Initialize controller setting PID parameters from kp, ki, and kd
    struct helm_state h;
    helm_reset(&h);
    h.kp = kp;         // Unified gain
    h.Td = kd / h.kp;  // Convert to derivative time scale
    h.Tf = f;          // Astrom and Murray p.308 suggests (h.Td / 2--20)
    h.Ti = h.kp / ki;  // Convert to integral time scale

    // Simulate controlled model, outputting status after each step
    // TODO Clean up the absolute time handling as roundoff avoidance got ugly
    helm_approach(&h);
    for (size_t i = 0; i*t < T+t; ++i) {
        v[0] += helm_steady(&h, t, r, u[0], v[0], y[0]);         // Control
        u[0]  = v[0];                                            // Ideal
        advance((i*t > T ? T - (i-1)*t : t), a, b, u, y);        // Advance
        printf("%-22.16g\t%-22.16g\t%-22.16g\t%-22.16g\t%-22.16g\n",
               i*t > T ? T : i*t, u[0], y[0], y[1], y[2]);       // Output
    }

    return EXIT_SUCCESS;
}
