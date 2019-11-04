/* This file is part of SWIFT.
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

/**
 *  @file mpipacked.c
 *  @brief file of routines to support packed particle exchanges.
 */

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <stddef.h>
#include <string.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local includes. */
#include "mpipacked.h"

#include "cell.h"
#include "part.h"

/* Whether the current "struct part" implementation supports packed
 * xv exchanges. */
int mpipacked_xvparts_supported = 0;

/* Whether the current "struct part" implementation supports packed
 * rho exchanges. */
int mpipacked_rhoparts_supported = 0;

/* Whether the current "struct gpart" implementation supports packed
 * exchanges. */
int mpipacked_gparts_supported = 0;

#ifdef WITH_MPI
/* The xv offsets and sizes, and number of these. */
static struct mpipacked_member *mpipacked_xvparts_members = NULL;
static int mpipacked_xvparts_nmembers = 0;

/**
 * @brief Create an MPI type to support the packed xv part exchanges.
 *
 * If mpipacked_xvparts_supported is false then this just returns a copy of
 * the full parts MPI type. Either way the new type should be committed
 * and freed.
 *
 * @param members array of member structs that defines the elements to extract.
 * @param nmembers number of elements in members array.
 * @param mpi_xvtype the new type.
 */
void mpipacked_make_type_xv(struct mpipacked_member members[],
                            int nmembers, MPI_Datatype *mpi_xvtype) {
  if (mpipacked_xvparts_supported) {
    size_t size = 0;
    for (int j = 0; j < nmembers; j++) {
      size += members[j].size;
    }
    if (MPI_Type_contiguous(size, MPI_BYTE, mpi_xvtype) != MPI_SUCCESS) {
      error("Failed to create a MPI type for xv updates");
    }

    /* Local references. */
    mpipacked_xvparts_members = members;
    mpipacked_xvparts_nmembers = nmembers;

  } else {

    /* So support for partial updates, so duplicate the existing type. */
    MPI_Type_dup(part_mpi_type, mpi_xvtype);
  }
}
#endif

/**
 * @brief Pack xv update members of parts in a cell.
 *
 * @param c The #cell.
 * @param packeddata an array to hold the bytes extracted from the parts
 *                   struct, must be large enough for all the extracted
 *                   members from all the part structs in the cell.
 */
void mpipacked_pack_parts_xv(struct cell *c, char *packeddata) {
#ifdef WITH_MPI
  size_t index = 0;
  for (int k = 0; k < c->hydro.count; k++) {
    char *part = (char *)&c->hydro.parts[k];
    for (int j = 0; j < mpipacked_xvparts_nmembers; j++) {
      memcpy(&packeddata[index], &part[mpipacked_xvparts_members[j].offset],
             mpipacked_xvparts_members[j].size);
      index += mpipacked_xvparts_members[j].size;
    }
  }
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Unpack the xv update members of parts to parts of a cell.
 *
 * @param c The #cell.
 * @param packeddata an array holding the bytes extracted from the parts
 *                   struct that need replacing.
 */
void mpipacked_unpack_parts_xv(struct cell *c, char *packeddata) {
#ifdef WITH_MPI
  size_t index = 0;
  for (int k = 0; k < c->hydro.count; k++) {
    char *part = (char *)&c->hydro.parts[k];
    for (int j = 0; j < mpipacked_xvparts_nmembers; j++) {
      memcpy(&part[mpipacked_xvparts_members[j].offset],
             &packeddata[index], mpipacked_xvparts_members[j].size);
      index += mpipacked_xvparts_members[j].size;
    }
  }
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
