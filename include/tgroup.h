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

#ifndef _tgroup_h
#define _tgroup_h
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"
#include "network.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
void init_ekindata(FILE *log,gmx_mtop_t *mtop,t_grpopts *opts,
			  gmx_ekindata_t *ekind);
=======
GMX_LIBMD_EXPORT
void init_ekindata(FILE *log, gmx_mtop_t *mtop, t_grpopts *opts,
                   gmx_ekindata_t *ekind);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Allocate memory and set the grpnr array. */

void done_ekindata(gmx_ekindata_t *ekind);
/* Free the memory */

<<<<<<< HEAD
void accumulate_u(t_commrec *cr,t_grpopts *opts,
			 gmx_ekindata_t *ekind);
=======
void accumulate_u(t_commrec *cr, t_grpopts *opts,
                  gmx_ekindata_t *ekind);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/*extern void accumulate_ekin(t_commrec *cr,t_grpopts *opts,t_groups *grps);*/
/* Communicate subsystem - group velocities and subsystem ekin respectively
 * and sum them up. Return them in grps.
 */

<<<<<<< HEAD
real sum_ekin(t_grpopts *opts,gmx_ekindata_t *ekind, real *dekindlambda, 
		     gmx_bool bEkinFullStep,gmx_bool bSaveEkinOld, gmx_bool bScaleEkin);
=======
GMX_LIBMD_EXPORT
real sum_ekin(t_grpopts *opts, gmx_ekindata_t *ekind, real *dekindlambda,
              gmx_bool bEkinFullStep, gmx_bool bSaveEkinOld, gmx_bool bScaleEkin);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Sum the group ekins into total ekin and calc temp per group,
 * return total temperature.
 */

<<<<<<< HEAD
void update_ekindata(int start,int homenr,gmx_ekindata_t *ekind,
			    t_grpopts *opts,rvec v[],t_mdatoms *md,real lambda);
=======
void update_ekindata(int start, int homenr, gmx_ekindata_t *ekind,
                     t_grpopts *opts, rvec v[], t_mdatoms *md, real lambda);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Do the update of group velocities (if bNEMD) and
 * (partial) group ekin.
 */

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _tgroup_h */
=======
#endif  /* _tgroup_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
