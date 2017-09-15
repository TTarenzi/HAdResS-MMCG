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

#ifndef _orires_h
#define _orires_h
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "sysstuff.h"
#include "typedefs.h"

#ifdef __cplusplus
<<<<<<< HEAD
    extern "C" {
#endif

void init_orires(FILE *fplog,const gmx_mtop_t *mtop,
			rvec x[],
			const t_inputrec *ir,
			const gmx_multisim_t *ms,t_oriresdata *od,
			t_state *state);
/* Initializes all the orientation restraint stuff in *od */

real calc_orires_dev(const gmx_multisim_t *ms,
			    int nfa,const t_iatom fa[],const t_iparams ip[],
			    const t_mdatoms *md,const rvec x[],
			    const t_pbc *pbc,t_fcdata *fcd,history_t *hist);
/* 
=======
extern "C" {
#endif

GMX_LIBGMX_EXPORT
void init_orires(FILE *fplog, const gmx_mtop_t *mtop,
                 rvec x[],
                 const t_inputrec *ir,
                 const gmx_multisim_t *ms, t_oriresdata *od,
                 t_state *state);
/* Initializes all the orientation restraint stuff in *od */

real calc_orires_dev(const gmx_multisim_t *ms,
                     int nfa, const t_iatom fa[], const t_iparams ip[],
                     const t_mdatoms *md, const rvec x[],
                     const t_pbc *pbc, t_fcdata *fcd, history_t *hist);
/*
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * Calculates the time averaged D matrices, the S matrix for each experiment.
 * Returns the weighted RMS deviation of the orientation restraints.
 */

<<<<<<< HEAD
=======
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
void diagonalize_orires_tensors(t_oriresdata *od);
/*
 * Diagonalizes the order tensor(s) of the orienation restraints.
 * For each experiment eig containts first 3 eigenvalues and then
 * the 3 eigenvectors. The eigenvalues are ordered on magnitude.
 */

<<<<<<< HEAD
void print_orires_log(FILE *log,t_oriresdata *od);
=======
GMX_LIBGMX_EXPORT
void print_orires_log(FILE *log, t_oriresdata *od);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Print order parameter, eigenvalues and eigenvectors to the log file */

t_ifunc orires;
/* Does only the orientation restraint force calculation */

<<<<<<< HEAD
void update_orires_history(t_fcdata *fcd,history_t *hist);
=======
GMX_LIBGMX_EXPORT
void update_orires_history(t_fcdata *fcd, history_t *hist);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Copy the new time averages that have been calculated in calc_orires_dev */

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _orires_h */
=======
#endif  /* _orires_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
