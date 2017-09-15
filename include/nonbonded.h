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

#ifndef _nonbonded_h
#define _nonbonded_h
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"
#include "pbc.h"
#include "network.h"
#include "tgroup.h"
#include "genborn.h"

#ifdef __cplusplus
extern "C" {
#endif
<<<<<<< HEAD

void gmx_setup_kernels(FILE *fplog,gmx_bool bGenericKernelOnly);
void gmx_setup_adress_kernels(FILE *fplog,gmx_bool bGenericKernelOnly);

#define GMX_DONB_LR             (1<<0)
#define GMX_DONB_FORCES         (1<<1)
#define GMX_DONB_FOREIGNLAMBDA  (1<<2)

void
do_nonbonded(t_commrec *cr,t_forcerec *fr,
             rvec x[],rvec f[],t_mdatoms *md,t_blocka *excl,
             real egnb[],real egcoul[],real egb[],rvec box_size,
             t_nrnb *nrnb,real lambda,real *dvdlambda,
             int nls,int eNL,int flags);

/* Calculate VdW/charge pair interactions (usually 1-4 interactions).
 * global_atom_index is only passed for printing error messages.
 */
real
do_listed_vdw_q(int ftype,int nbonds,
		const t_iatom iatoms[],const t_iparams iparams[],
		const rvec x[],rvec f[],rvec fshift[],
		const t_pbc *pbc,const t_graph *g,
		real lambda,real *dvdlambda,
		const t_mdatoms *md,
		const t_forcerec *fr,gmx_grppairener_t *grppener,
		int *global_atom_index);


void
adress_drift_term(FILE *               fplog,
                          int                  cg0,
                          int                  cg1,
                          int                   cg1home,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc, rvec f[], real *engdelta);


=======
#if 0
} /* fixes auto-indentation problems */
#endif



GMX_LIBGMX_EXPORT
void
gmx_nonbonded_setup(FILE *         fplog,
                    t_forcerec *   fr,
                    gmx_bool       bGenericKernelOnly);





GMX_LIBGMX_EXPORT
void
gmx_nonbonded_set_kernel_pointers(FILE *       fplog,
                                  t_nblist *   nl);



#define GMX_NONBONDED_DO_LR             (1<<0)
#define GMX_NONBONDED_DO_FORCE          (1<<1)
#define GMX_NONBONDED_DO_FOREIGNLAMBDA  (1<<2)
#define GMX_NONBONDED_DO_POTENTIAL      (1<<3)
#define GMX_NONBONDED_DO_SR             (1<<4)

GMX_LIBGMX_EXPORT
void
do_nonbonded(t_commrec *cr, t_forcerec *fr,
             rvec x[], rvec f_shortrange[], rvec f_longrange[], t_mdatoms *md, t_blocka *excl,
             gmx_grppairener_t *grppener, rvec box_size,
             t_nrnb *nrnb, real *lambda, real dvdlambda[],
             int nls, int eNL, int flags);

/* Calculate VdW/charge listed pair interactions (usually 1-4 interactions).
 * global_atom_index is only passed for printing error messages.
 */
real
do_nonbonded_listed(int ftype, int nbonds, const t_iatom iatoms[], const t_iparams iparams[],
                    const rvec x[], rvec f[], rvec fshift[], const t_pbc *pbc, const t_graph *g,
                    real *lambda, real *dvdl, const t_mdatoms *md, const t_forcerec *fr,
                    gmx_grppairener_t *grppener, int *global_atom_index);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef __cplusplus
}
#endif

#endif
