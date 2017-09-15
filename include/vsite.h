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

#ifndef _vsite_h
#define _vsite_h
<<<<<<< HEAD

#include <stdio.h>
#include "typedefs.h"
=======
#include <stdio.h>
#include "visibility.h"
#include "typedefs.h"
#include "types/commrec.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
<<<<<<< HEAD
        int *      left_import_construct;
        int        left_import_nconstruct;
        int *      left_export_construct;
        int        left_export_nconstruct;
        int *      right_import_construct;
        int        right_import_nconstruct;
        int *      right_export_construct;
        int        right_export_nconstruct;
        rvec *     send_buf;
        rvec *     recv_buf;
} t_comm_vsites;

typedef struct {
  int  n_intercg_vsite;       /* The number of inter charge group vsites */
  int  nvsite_pbc_molt;       /* The array size of vsite_pbc_molt        */
  int  ***vsite_pbc_molt;     /* The pbc atoms for intercg vsites        */
  int  **vsite_pbc_loc;       /* The local pbc atoms                     */
  int  *vsite_pbc_loc_nalloc;
  gmx_bool bPDvsitecomm;          /* Do we need vsite communication with PD? */
  t_comm_vsites *vsitecomm;   /* The PD vsite communication struct       */
} gmx_vsite_t;

void construct_vsites(FILE *log,gmx_vsite_t *vsite,
			     rvec x[],t_nrnb *nrnb,
			     real dt,rvec v[],
			     t_iparams ip[],t_ilist ilist[],
			     int ePBC,gmx_bool bMolPBC,t_graph *graph,
			     t_commrec *cr,matrix box);
=======
    int *      left_import_construct;
    int        left_import_nconstruct;
    int *      left_export_construct;
    int        left_export_nconstruct;
    int *      right_import_construct;
    int        right_import_nconstruct;
    int *      right_export_construct;
    int        right_export_nconstruct;
    rvec *     send_buf;
    rvec *     recv_buf;
} t_comm_vsites;

typedef struct {
    t_ilist ilist[F_NRE];     /* vsite ilists for this thread            */
    rvec    fshift[SHIFTS];   /* fshift accumulation buffer              */
    matrix  dxdf;             /* virial dx*df accumulation buffer        */
} gmx_vsite_thread_t;

typedef struct {
    gmx_bool            bHaveChargeGroups;    /* Do we have charge groups?               */
    int                 n_intercg_vsite;      /* The number of inter charge group vsites */
    int                 nvsite_pbc_molt;      /* The array size of vsite_pbc_molt        */
    int              ***vsite_pbc_molt;       /* The pbc atoms for intercg vsites        */
    int               **vsite_pbc_loc;        /* The local pbc atoms                     */
    int                *vsite_pbc_loc_nalloc; /* Sizes of vsite_pbc_loc                  */
    gmx_bool            bPDvsitecomm;         /* Do we need vsite communication with PD? */
    t_comm_vsites      *vsitecomm;            /* The PD vsite communication struct       */
    int                 nthreads;             /* Number of threads used for vsites       */
    gmx_vsite_thread_t *tdata;                /* Thread local vsites and work structs    */
    int                *th_ind;               /* Work array                              */
    int                 th_ind_nalloc;        /* Size of th_ind                          */
} gmx_vsite_t;

GMX_LIBMD_EXPORT
void construct_vsites(FILE *log, gmx_vsite_t *vsite,
                      rvec x[], t_nrnb *nrnb,
                      real dt, rvec v[],
                      t_iparams ip[], t_ilist ilist[],
                      int ePBC, gmx_bool bMolPBC, t_graph *graph,
                      t_commrec *cr, matrix box);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Create positions of vsite atoms based on surrounding atoms
 * for the local system.
 * If v is passed, the velocities of the vsites will be calculated
 * as the new positions minus the old positions divided by dt,
 * thus v should only be passed when the coordinates have been
 * updated with a full time step.
 * Note that velocitis of vsites are completely irrelevant
 * for the integration, they are only useful for analysis.
 */

<<<<<<< HEAD
void construct_vsites_mtop(FILE *log,gmx_vsite_t *vsite,
			   gmx_mtop_t *mtop,rvec x[]);
=======
GMX_LIBMD_EXPORT
void construct_vsites_mtop(FILE *log, gmx_vsite_t *vsite,
                           gmx_mtop_t *mtop, rvec x[]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Create positions of vsite atoms based on surrounding atoms
 * for the whole system.
 * This function assumes that all molecules are whole.
 */

<<<<<<< HEAD
void spread_vsite_f(FILE *log,gmx_vsite_t *vsite,
			   rvec x[],rvec f[],rvec *fshift,
			   t_nrnb *nrnb,t_idef *idef,
			   int ePBC,gmx_bool bMolPBC,t_graph *g,matrix box,
			   t_commrec *cr);
/* Spread the force operating on the vsite atoms on the surrounding atoms.
 * If fshift!=NULL also update the shift forces.
 */

gmx_vsite_t *init_vsite(gmx_mtop_t *mtop,t_commrec *cr);
/* Initialize the virtual site struct,
 * returns NULL when there are no virtual sites.
 */

void set_vsite_top(gmx_vsite_t *vsite,gmx_localtop_t *top,t_mdatoms *md,
			  t_commrec *cr);
=======
void spread_vsite_f(FILE *log, gmx_vsite_t *vsite,
                    rvec x[], rvec f[], rvec *fshift,
                    gmx_bool VirCorr, matrix vir,
                    t_nrnb *nrnb, t_idef *idef,
                    int ePBC, gmx_bool bMolPBC, t_graph *g, matrix box,
                    t_commrec *cr);
/* Spread the force operating on the vsite atoms on the surrounding atoms.
 * If fshift!=NULL also update the shift forces.
 * If VirCorr=TRUE add the virial correction for non-linear vsite constructs
 * to vir. This correction is required when the virial is not calculated
 * afterwards from the particle position and forces, but in a different way,
 * as for instance for the PME mesh contribution.
 */

GMX_LIBMD_EXPORT
gmx_vsite_t *init_vsite(gmx_mtop_t *mtop, t_commrec *cr,
                        gmx_bool bSerial_NoPBC);
/* Initialize the virtual site struct,
 * returns NULL when there are no virtual sites.
 * bSerial_NoPBC is to generate a simple vsite setup to be
 * used only serial (no MPI or thread parallelization) and without pbc;
 * this is useful for correction vsites of the initial configuration.
 */

void split_vsites_over_threads(const t_ilist   *ilist,
                               const t_mdatoms *mdatoms,
                               gmx_bool         bLimitRange,
                               gmx_vsite_t     *vsite);
/* Divide the vsite work-load over the threads.
 * Should be called at the end of the domain decomposition.
 */

GMX_LIBMD_EXPORT
void set_vsite_top(gmx_vsite_t *vsite, gmx_localtop_t *top, t_mdatoms *md,
                   t_commrec *cr);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Set some vsite data for runs without domain decomposition.
 * Should be called once after init_vsite, before calling other routines.
 */

#ifdef __cplusplus
}
#endif

#endif
<<<<<<< HEAD

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2