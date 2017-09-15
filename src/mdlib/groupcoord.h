<<<<<<< HEAD
/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */

/*! \file groupcoord.h
 *
 *  @brief Assemble atom positions for comparison with a reference set.
 *
<<<<<<< HEAD
 *  This file contains functions to assemble the positions of a subset of the 
=======
 *  This file contains functions to assemble the positions of a subset of the
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *  atoms and to do operations on it like determining the center of mass, or
 *  doing translations and rotations. These functions are useful when
 *  a subset of the positions needs to be compared to some set of reference
 *  positions, as e.g. done for essential dynamics.
<<<<<<< HEAD
 *  
=======
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
<<<<<<< HEAD


/*! \brief Select local atoms of a group.
*
* Selects the indices of local atoms of a group and stores them in anrs_loc[0..nr_loc]. 
* If you need the positions of the group's atoms on all nodes, provide a coll_ind[0..nr] 
* array and pass it on to communicate_group_positions. Thus the collective array 
* will always have the same atom order (ascending indices). 
*
*  \param ga2la[in]         Global to local atom index conversion data.
*  \param nr[in]            The total number of atoms that the group contains.
*  \param anrs[in]          The global atom number of the group's atoms.
*  \param nr_loc[out]       The number of group atoms present on the local node.
*  \param anrs_loc[out]     The local atom numbers of the group.
*  \param nalloc_loc[inout] Local allocation size of anrs_loc array.
*  \param coll_ind[opt]     If not NULL this array must be of size nr. It stores
*                           for each local atom where it belongs in the global
*                           (collective) array such that it can be gmx_summed
*                           in the communicate_group_positions routine.
*/
=======
#include "types/commrec.h"


/*! \brief Select local atoms of a group.
 *
 * Selects the indices of local atoms of a group and stores them in anrs_loc[0..nr_loc].
 * If you need the positions of the group's atoms on all nodes, provide a coll_ind[0..nr]
 * array and pass it on to communicate_group_positions. Thus the collective array
 * will always have the same atom order (ascending indices).
 *
 *  \param ga2la[in]         Global to local atom index conversion data.
 *  \param nr[in]            The total number of atoms that the group contains.
 *  \param anrs[in]          The global atom number of the group's atoms.
 *  \param nr_loc[out]       The number of group atoms present on the local node.
 *  \param anrs_loc[out]     The local atom numbers of the group.
 *  \param nalloc_loc[inout] Local allocation size of anrs_loc array.
 *  \param coll_ind[opt]     If not NULL this array must be of size nr. It stores
 *                           for each local atom where it belongs in the global
 *                           (collective) array such that it can be gmx_summed
 *                           in the communicate_group_positions routine.
 */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
extern void dd_make_local_group_indices(gmx_ga2la_t ga2la,
                                        const int nr, int anrs[], int *nr_loc,
                                        int *anrs_loc[], int *nalloc_loc,
                                        int coll_ind[]);


/*! \brief Assemble local positions into a collective array present on all nodes.
<<<<<<< HEAD
 * 
 * Communicate the positions of the group's atoms such that every node has all of 
 * them. Unless running on huge number of cores, this is not a big performance impact
 * as long as the collective subset [0..nr] is kept small. The atom indices are 
 * retrieved from anrs_loc[0..nr_loc]. If you call the routine for the serial case,
 * provide an array coll_ind[i] = i for i in 1..nr.
 * 
 * \param cr[in]             Pointer to MPI communication data.     
=======
 *
 * Communicate the positions of the group's atoms such that every node has all of
 * them. Unless running on huge number of cores, this is not a big performance impact
 * as long as the collective subset [0..nr] is kept small. The atom indices are
 * retrieved from anrs_loc[0..nr_loc]. If you call the routine for the serial case,
 * provide an array coll_ind[i] = i for i in 1..nr.
 *
 * \param cr[in]             Pointer to MPI communication data.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * \param xcoll[out]         Collective array of positions, idential on all nodes
 *                           after this routine has been called.
 * \param shifts[inout]      Collective array of shifts for xcoll, needed to make
 *                           the group whole. This array remembers the shifts
 *                           since the start of the simulation (where the group
 *                           is whole) and must therefore not be changed outside
 *                           of this routine!
 * \param extra_shifts[buf]  Extra shifts since last time step, only needed as
 *                           buffer variable [0..nr].
<<<<<<< HEAD
 * \param bNS[in]            Neighborsearching/domain redecomposition has been 
 *                           performed at the begin of this time step such that 
 *                           the shifts have changed and need to be updated.
 * \param x_loc[in]          Pointer to the local atom positions this node has. 
 * \param nr[in]             Total number of atoms in the group.
 * \param nr_loc[in]         Number of group atoms on the local node.
 * \param anrs_loc[in]       Array of the local atom indices.
 * \param coll_ind[in]       This array of size nr stores for each local atom where 
=======
 * \param bNS[in]            Neighborsearching/domain redecomposition has been
 *                           performed at the begin of this time step such that
 *                           the shifts have changed and need to be updated.
 * \param x_loc[in]          Pointer to the local atom positions this node has.
 * \param nr[in]             Total number of atoms in the group.
 * \param nr_loc[in]         Number of group atoms on the local node.
 * \param anrs_loc[in]       Array of the local atom indices.
 * \param coll_ind[in]       This array of size nr stores for each local atom where
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *                           it belongs in the collective array so that the local
 *                           contributions can be gmx_summed. It is provided by
 *                           dd_make_local_group_indices.
 * \param xcoll_old[inout]   Positions from the last time step, used to make the
 *                           group whole.
 * \param box[in]            Simulation box matrix, needed to shift xcoll such that
 *                           the group becomes whole.
 */
extern void communicate_group_positions(t_commrec *cr, rvec *xcoll, ivec *shifts,
                                        ivec *extra_shifts, const gmx_bool bNS,
                                        rvec *x_loc, const int nr, const int nr_loc,
                                        int *anrs_loc, int *coll_ind, rvec *xcoll_old,
                                        matrix box);


/*! \brief Calculates the center of the positions x locally.
<<<<<<< HEAD
 * 
 * Calculates the center of mass (if masses are given in the weight array) or
 * the geometrical center (if NULL is passed as weight).
 * 
=======
 *
 * Calculates the center of mass (if masses are given in the weight array) or
 * the geometrical center (if NULL is passed as weight).
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * \param x[in]              Positions.
 * \param weight[in]         Can be NULL or an array of weights. If masses are
 *                           given as weights, the COM is calculated.
 * \param nr[in]             Number of positions and weights if present.
 * \param center[out]        The (weighted) center of the positions.
<<<<<<< HEAD
 * 
=======
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */
extern void get_center(rvec x[], real weight[], const int nr, rvec center);


/*! \brief Calculates the sum of the positions x locally.
<<<<<<< HEAD
 * 
 * Calculates the (weighted) sum of position vectors and returns the sum of 
 * weights, which is needed when local contributions shall be summed to a 
 * global weighted center.
 * 
=======
 *
 * Calculates the (weighted) sum of position vectors and returns the sum of
 * weights, which is needed when local contributions shall be summed to a
 * global weighted center.
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * \param x[in]              Array of positions.
 * \param weight[in]         Can be NULL or an array of weights.
 * \param nr[in]             Number of positions and weights if present.
 * \param dsumvec[out]       The (weighted) sum of the positions.
 * \return Sum of weights.
<<<<<<< HEAD
 * 
=======
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */
extern double get_sum_of_positions(rvec x[], real weight[], const int nr, dvec dsumvec);


<<<<<<< HEAD
/*! \brief Calculates the global center of all local arrays x. 
 * 
 * Get the center from local positions [0..nr_loc], this involves communication.
 * Not that the positions must already have the correct PBC representation. Use
 * this routine if no collective coordinates are assembled from which the center 
 * could be calculated without communication.
 * 
=======
/*! \brief Calculates the global center of all local arrays x.
 *
 * Get the center from local positions [0..nr_loc], this involves communication.
 * Not that the positions must already have the correct PBC representation. Use
 * this routine if no collective coordinates are assembled from which the center
 * could be calculated without communication.
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * \param cr[in]             Pointer to MPI communication data.
 * \param x_loc[in]          Array of local positions [0..nr_loc].
 * \param weight_loc[in]     Array of local weights, these are the masses if the
 *                           center of mass is to be calculated.
 * \param nr_loc[in]         The number of positions on the local node.
<<<<<<< HEAD
 * \param nr_group[in]       The number of positions in the whole group. Since 
=======
 * \param nr_group[in]       The number of positions in the whole group. Since
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *                           this is known anyway, we do not need to communicate
 *                           and sum nr_loc if we pass it over.
 * \param center[out]        The (weighted) center of all x_loc from all the
 *                           nodes.
 */
extern void get_center_comm(t_commrec *cr, rvec x_loc[], real weight_loc[],
                            int nr_loc, int nr_group, rvec center);


/*! \brief Translate positions.
<<<<<<< HEAD
 * 
 * Add a translation vector to the positions x.
 * 
 * \param x[inout]           Array of positions.
 * \param nr[in]             Number of entries in the position array.
 * \param transvec[in]       Translation vector to be added to all positions.
 * 
=======
 *
 * Add a translation vector to the positions x.
 *
 * \param x[inout]           Array of positions.
 * \param nr[in]             Number of entries in the position array.
 * \param transvec[in]       Translation vector to be added to all positions.
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */
extern void translate_x(rvec x[], const int nr, const rvec transvec);


/*! \brief Rotate positions.
<<<<<<< HEAD
 * 
 * Rotate the positions with the rotation matrix.
 * 
 * \param x[inout]           Array of positions.
 * \param nr[in]             Number of entries in the position array.
 * \param rmat[in]           Rotation matrix to operate on all positions.
 * 
 */
extern void rotate_x(rvec x[], const int nr, matrix rmat);

=======
 *
 * Rotate the positions with the rotation matrix.
 *
 * \param x[inout]           Array of positions.
 * \param nr[in]             Number of entries in the position array.
 * \param rmat[in]           Rotation matrix to operate on all positions.
 *
 */
extern void rotate_x(rvec x[], const int nr, matrix rmat);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
