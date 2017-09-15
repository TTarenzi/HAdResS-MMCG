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

#ifndef _mvdata_h
#define _mvdata_h
<<<<<<< HEAD

#include "typedefs.h"
=======
#include "visibility.h"
#include "typedefs.h"
#include "types/commrec.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
void bcast_ir_mtop(const t_commrec *cr,
			 t_inputrec *inputrec,gmx_mtop_t *mtop);
/* Broadcasts ir and mtop from the master to all nodes in cr->mpi_comm_mygroup.
 */

void bcast_state_setup(const t_commrec *cr,t_state *state);
=======
GMX_LIBGMX_EXPORT
void bcast_ir_mtop(const t_commrec *cr,
                   t_inputrec *inputrec, gmx_mtop_t *mtop);
/* Broadcasts ir and mtop from the master to all nodes in cr->mpi_comm_mygroup.
 */

GMX_LIBGMX_EXPORT
void bcast_state_setup(const t_commrec *cr, t_state *state);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Broadcasts the state sizes and flags
 * from the master to all nodes in cr->mpi_comm_mygroup.
 * The arrays are not broadcasted.
 */

<<<<<<< HEAD
void bcast_state(const t_commrec *cr,t_state *state,gmx_bool bAlloc);
=======
GMX_LIBGMX_EXPORT
void bcast_state(const t_commrec *cr, t_state *state, gmx_bool bAlloc);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Broadcasts state from the master to all nodes in cr->mpi_comm_mygroup.
 * The arrays in state are allocated when bAlloc is TRUE.
 */


/* Routines for particle decomposition only in mvxvf.c */

<<<<<<< HEAD
void move_cgcm(FILE *log,const t_commrec *cr,rvec cg_cm[]);
		     
void move_rvecs(const t_commrec *cr,gmx_bool bForward,gmx_bool bSum,
		       int left,int right,rvec vecs[],rvec buf[],
		       int shift,t_nrnb *nrnb);

void move_reals(const t_commrec *cr,gmx_bool bForward,gmx_bool bSum,
                       int left,int right,real reals[],real buf[],
                       int shift,t_nrnb *nrnb);

void move_x(FILE *log,const t_commrec *cr,
		   int left,int right,rvec x[],t_nrnb *nrnb);
		    
void move_rborn(FILE *log,const t_commrec *cr,
                       int left,int right,real rborn[],t_nrnb *nrnb);

void move_f(FILE *log,const t_commrec *cr,
		   int left,int right,rvec f[],rvec fadd[],
		   t_nrnb *nrnb);

void move_gpol(FILE *log,const t_commrec *cr,
                      int left,int right,real gpol[],real gpol_add[],
                      t_nrnb *nrnb);
=======
void move_cgcm(FILE *log, const t_commrec *cr, rvec cg_cm[]);

GMX_LIBMD_EXPORT
void move_rvecs(const t_commrec *cr, gmx_bool bForward, gmx_bool bSum,
                int left, int right, rvec vecs[], rvec buf[],
                int shift, t_nrnb *nrnb);

void move_reals(const t_commrec *cr, gmx_bool bForward, gmx_bool bSum,
                int left, int right, real reals[], real buf[],
                int shift, t_nrnb *nrnb);

void move_x(FILE *log, const t_commrec *cr,
            int left, int right, rvec x[], t_nrnb *nrnb);

void move_rborn(FILE *log, const t_commrec *cr,
                int left, int right, real rborn[], t_nrnb *nrnb);

void move_f(FILE *log, const t_commrec *cr,
            int left, int right, rvec f[], rvec fadd[],
            t_nrnb *nrnb);

void move_gpol(FILE *log, const t_commrec *cr,
               int left, int right, real gpol[], real gpol_add[],
               t_nrnb *nrnb);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _mvdata_h */
=======
#endif  /* _mvdata_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
