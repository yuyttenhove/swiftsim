/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_RT_ADDITIONS_H
#define SWIFT_RT_ADDITIONS_H

/**
 * @file src/rt_additions.h
 * @brief Branches between the different additional functions required outside
 * of the RT files; Required to avoid circular inclusions.
 */

/* Config parameters. */
#include "../config.h"

/* Import the right RT definition */
#if defined(RT_NONE)
#include "./rt/none/rt_additions.h"
#elif defined(RT_DEBUG)
#include "./rt/debug/rt_additions.h"
#elif defined(RT_GEAR)
#include "./rt/GEAR/rt_additions.h"
#elif defined(RT_SPHM1RT)
#include "./rt/SPHM1RT/rt_additions.h"
#else
#error "Invalid choice of radiation scheme"
#endif

#endif /* SWIFT_RT_ADDITIONS_H */
