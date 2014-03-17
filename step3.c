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
 * Example controlling a 3rd order ODE across a unit step in the setpoint.
 *
 * FIXME
 * This system of ODEs may be used to simulate the temporal response of a system
 * given a unit step input.  That is, just prior to time zero process state
 * \f$y(t)\f$, reference value \f$r(t)\f$, actuator signal \f$u(t)\f$, and all
 * of their derivatives are zero.  At time zero, step change \f$r(t) = 1\f$ is
 * introduced.  These ODEs, in conjunction with the controller dynamics (and
 * any discretization choices), determine the controlled system response.
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "helm.h"

/**
 * Advance the temporal state of a model given by transfer function
 * \f$ \frac{y(s)}{u(s)} = \frac{b_0}{s^3 + a_2 s^2 + a_1 s + a_0} \f$.
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
        const double a[static 3],
        const double b[static 1],
        const double u[static 1],
              double y[static 3]);

int
main (int argc, char *argv[])
{
    // TODO
    double a[3] = {1, 3, 3};  // Based on Astrom & Murray Figure 10.2 example
    double b[1] = {1};        // Ditto

    return EXIT_SUCCESS;
}
