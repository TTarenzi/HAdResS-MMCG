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

#ifndef _ns_h
#define _ns_h

#include <stdio.h>
<<<<<<< HEAD
=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "sysstuff.h"
#include "typedefs.h"
#include "pbc.h"
#include "tgroup.h"
#include "network.h"


#ifdef __cplusplus
extern "C" {
#endif

/****************************************************
 *
 *    U T I L I T I E S May be found in ns.c
 *
 ****************************************************/

<<<<<<< HEAD
void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr);
/* 
=======
GMX_LIBMD_EXPORT
void init_neighbor_list(FILE *log, t_forcerec *fr, int homenr);
/*
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * nn is the number of energy terms in the energy matrix
 * (ngener*(ngener-1))/2
 * start is the first atom on this processor
 * homenr is the number of atoms on this processor
 */
<<<<<<< HEAD
 
int calc_naaj(int icg,int cgtot);
=======

int calc_naaj(int icg, int cgtot);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Calculate the number of charge groups to interact with for icg */

/****************************************************
 *
 *    N E I G H B O R  S E A R C H I N G
 *
 *    Calls either ns5_core (when grid selected in .mdp file)
 *    or ns_simple_core (when simple selected in .mdp file)
 *
<<<<<<< HEAD
 *    Return total number of pairs searched 
 *
 ****************************************************/
void init_ns(FILE *fplog,const t_commrec *cr,
		    gmx_ns_t *ns,t_forcerec *fr,
		    const gmx_mtop_t *mtop,
		    matrix box);

int search_neighbours(FILE *log,t_forcerec *fr,
			     rvec x[],matrix box,
			     gmx_localtop_t *top,
			     gmx_groups_t *groups,
			     t_commrec *cr,
			     t_nrnb *nrnb,t_mdatoms *md,
			     real lambda,real *dvdlambda,
			     gmx_grppairener_t *grppener,
			     gmx_bool bFillGrid,
			     gmx_bool bDoLongRange,
			     gmx_bool bDoForces,rvec *f);
 

/* Debugging routines from wnblist.c */
void dump_nblist(FILE *out,t_commrec *cr,t_forcerec *fr,int nDNL);

int read_nblist(FILE *in,FILE *out,int **mat,int natoms,gmx_bool bSymm);
/* Returns total number of neighbors. If bSymm the matrix is symmetrized. */

int natoms_beyond_ns_buffer(t_inputrec *ir,t_forcerec *fr,t_block *cgs,
				   matrix scale_tot,rvec *x);
=======
 *    Return total number of pairs searched
 *
 ****************************************************/
void init_ns(FILE *fplog, const t_commrec *cr,
             gmx_ns_t *ns, t_forcerec *fr,
             const gmx_mtop_t *mtop,
             matrix box);

GMX_LIBMD_EXPORT
int search_neighbours(FILE *log, t_forcerec *fr,
                      rvec x[], matrix box,
                      gmx_localtop_t *top,
                      gmx_groups_t *groups,
                      t_commrec *cr,
                      t_nrnb *nrnb, t_mdatoms *md,
                      real *lambda, real *dvdlambda,
                      gmx_grppairener_t *grppener,
                      gmx_bool bFillGrid,
                      gmx_bool bDoLongRangeNS,
                      gmx_bool bPadListsForKernels);


/* Debugging routines from wnblist.c */
GMX_LIBMD_EXPORT
void dump_nblist(FILE *out, t_commrec *cr, t_forcerec *fr, int nDNL);

int read_nblist(FILE *in, FILE *out, int **mat, int natoms, gmx_bool bSymm);
/* Returns total number of neighbors. If bSymm the matrix is symmetrized. */

GMX_LIBMD_EXPORT
int natoms_beyond_ns_buffer(t_inputrec *ir, t_forcerec *fr, t_block *cgs,
                            matrix scale_tot, rvec *x);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Returns the number of atoms that moved beyond the ns buffer */

#ifdef __cplusplus
}
#endif


<<<<<<< HEAD
#endif	/* _ns_h */
=======
#endif  /* _ns_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
