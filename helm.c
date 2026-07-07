//--------------------------------------------------------------------------
//
// Copyright (C) 2014, 2026 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** \file
 * C99 extern declarations for static inline functions within \ref helm.h
 *
 * When a header-only library uses static inline, each translation unit
 * gets its own private copy.  Providing matching extern declarations in
 * a single .c file lets the compiler emit one shared, externally-visible
 * definition that the linker can use when inlining is not performed.
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
