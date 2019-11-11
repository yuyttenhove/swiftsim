/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_MPIPACKED_H
#define SWIFT_MPIPACKED_H

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <stddef.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Forward declarations. */
struct cell;

/* Whether the implementation supports packed exchanges. Bit mask with
 * task_mpitype_xxx values. */
extern int mpipacked_supported;

/* The byte offset and size of a struct member. We need a list of these to
 * support each packed exchange type. */
struct mpipacked_member {
  size_t offset;
  size_t size;
};

/* Macro that expands into a single members struct. */
#define MPIPACKED_ADDMEMBER(mystruct, member) \
  { offsetof(mystruct, member), sizeof(((mystruct *)0)->member) }

#ifdef WITH_MPI
void mpipacked_make_type_packed_xv(MPI_Datatype *mpi_type);
void mpipacked_make_type_packed_gpart(MPI_Datatype *mpi_type);
#endif
void mpipacked_pack_parts_xv(struct cell *c, char *packeddata);
void mpipacked_unpack_parts_xv(struct cell *c, char *packeddata);

void mpipacked_pack_gparts(struct cell *c, char *packeddata);
void mpipacked_unpack_gparts(struct cell *c, char *packeddata);

#endif /* SWIFT_MPIPACKED_H */
