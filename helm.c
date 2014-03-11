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
 * C99 extern declarations for state inline functions within \ref helm.h
 */

#include "helm.h"

extern
struct helm_state *
helm_reset(struct helm_state * const h);

extern
struct helm_state *
helm_approach(struct helm_state * const h);

extern
double
helm_steady(struct helm_state * const h,
            const double dt,
            const double r,
            const double u,
            const double v,
            const double y);
