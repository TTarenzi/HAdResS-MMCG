/*
<<<<<<< HEAD
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gromacs Runs On Most of All Computer Systems
=======
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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

#ifndef _princ_h
#define _princ_h
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
void rotate_atoms(int gnx,atom_id index[],rvec x[],matrix trans);
/* Rotate all atoms in index using matrix trans */

void principal_comp(int n,atom_id index[],t_atom atom[],rvec x[],
			   matrix trans,rvec d);
=======
void rotate_atoms(int gnx, atom_id index[], rvec x[], matrix trans);
/* Rotate all atoms in index using matrix trans */

GMX_LIBGMX_EXPORT
void principal_comp(int n, atom_id index[], t_atom atom[], rvec x[],
                    matrix trans, rvec d);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Calculate the principal components of atoms in index. Atoms are
 * mass weighted. It is assumed that the center of mass is in the origin!
 */

<<<<<<< HEAD
void orient_princ(t_atoms *atoms, int isize, atom_id *index,
			 int natoms, rvec x[], rvec *v, rvec d);
/* rotates molecule to align principal axes with coordinate axes */

real calc_xcm(rvec x[],int gnx,atom_id *index,t_atom *atom,rvec xcm,
		     gmx_bool bQ);
=======
GMX_LIBGMX_EXPORT
void orient_princ(t_atoms *atoms, int isize, atom_id *index,
                  int natoms, rvec x[], rvec *v, rvec d);
/* rotates molecule to align principal axes with coordinate axes */

GMX_LIBGMX_EXPORT
real calc_xcm(rvec x[], int gnx, atom_id *index, t_atom *atom, rvec xcm,
              gmx_bool bQ);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Calculate the center of mass of the atoms in index. if bQ then the atoms
 * will be charge weighted rather than mass weighted.
 * Returns the total mass/charge.
 */

<<<<<<< HEAD
real sub_xcm(rvec x[],int gnx,atom_id *index,t_atom atom[],rvec xcm,
		    gmx_bool bQ);
=======
GMX_LIBGMX_EXPORT
real sub_xcm(rvec x[], int gnx, atom_id *index, t_atom atom[], rvec xcm,
             gmx_bool bQ);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Calc. the center of mass and subtract it from all coordinates.
 * Returns the original center of mass in xcm
 * Returns the total mass
 */

<<<<<<< HEAD
void add_xcm(rvec x[],int gnx,atom_id *index,rvec xcm);
=======
void add_xcm(rvec x[], int gnx, atom_id *index, rvec xcm);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Increment all atoms in index with xcm */

#ifdef __cplusplus
}
#endif

#endif
