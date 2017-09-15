/*
<<<<<<< HEAD
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "domdec.h"
#include "gmx_wallcycle.h"
<<<<<<< HEAD
=======
#include "gmx_cyclecounter.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "trnio.h"
#include "smalloc.h"
#include "network.h"
#include "pbc.h"
#include "futil.h"
#include "mdrun.h"
#include "txtdump.h"
#include "names.h"
#include "mtop_util.h"
#include "names.h"
#include "nrjac.h"
#include "vec.h"
#include "gmx_ga2la.h"
#include "xvgr.h"
#include "gmxfio.h"
#include "groupcoord.h"
#include "pull_rotation.h"
#include "gmx_sort.h"
#include "copyrite.h"
<<<<<<< HEAD
#include "gmx_cyclecounter.h"
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


static char *RotStr = {"Enforced rotation:"};


/* Set the minimum weight for the determination of the slab centers */
#define WEIGHT_MIN (10*GMX_FLOAT_MIN)

/* Helper structure for sorting positions along rotation vector             */
typedef struct {
    real xcproj;            /* Projection of xc on the rotation vector        */
<<<<<<< HEAD
    int ind;                /* Index of xc                                    */
=======
    int  ind;               /* Index of xc                                    */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    real m;                 /* Mass                                           */
    rvec x;                 /* Position                                       */
    rvec x_ref;             /* Reference position                             */
} sort_along_vec_t;


/* Enforced rotation / flexible: determine the angle of each slab             */
typedef struct gmx_slabdata
{
<<<<<<< HEAD
    int  nat;               /* Number of atoms belonging to this slab         */
    rvec *x;                /* The positions belonging to this slab. In 
                               general, this should be all positions of the 
                               whole rotation group, but we leave those away 
=======
    int   nat;              /* Number of atoms belonging to this slab         */
    rvec *x;                /* The positions belonging to this slab. In
                               general, this should be all positions of the
                               whole rotation group, but we leave those away
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                               that have a small enough weight                */
    rvec *ref;              /* Same for reference                             */
    real *weight;           /* The weight for each atom                       */
} t_gmx_slabdata;


/* Helper structure for potential fitting */
typedef struct gmx_potfit
{
    real   *degangle;       /* Set of angles for which the potential is
                               calculated. The optimum fit is determined as
                               the angle for with the potential is minimal    */
    real   *V;              /* Potential for the different angles             */
    matrix *rotmat;         /* Rotation matrix corresponding to the angles    */
} t_gmx_potfit;


/* Enforced rotation data for all groups                                      */
typedef struct gmx_enfrot
{
<<<<<<< HEAD
    FILE  *out_rot;         /* Output file for rotation data                  */
    FILE  *out_torque;      /* Output file for torque data                    */
    FILE  *out_angles;      /* Output file for slab angles for flexible type  */
    FILE  *out_slabs;       /* Output file for slab centers                   */
    int   bufsize;          /* Allocation size of buf                         */
    rvec  *xbuf;            /* Coordinate buffer variable for sorting         */
    real  *mbuf;            /* Masses buffer variable for sorting             */
    sort_along_vec_t *data; /* Buffer variable needed for position sorting    */
    real  *mpi_inbuf;       /* MPI buffer                                     */
    real  *mpi_outbuf;      /* MPI buffer                                     */
    int   mpi_bufsize;      /* Allocation size of in & outbuf                 */
    unsigned long Flags;    /* mdrun flags                                    */
    gmx_bool bOut;          /* Used to skip first output when appending to 
                             * avoid duplicate entries in rotation outfiles   */
=======
    FILE             *out_rot;     /* Output file for rotation data                  */
    FILE             *out_torque;  /* Output file for torque data                    */
    FILE             *out_angles;  /* Output file for slab angles for flexible type  */
    FILE             *out_slabs;   /* Output file for slab centers                   */
    int               bufsize;     /* Allocation size of buf                         */
    rvec             *xbuf;        /* Coordinate buffer variable for sorting         */
    real             *mbuf;        /* Masses buffer variable for sorting             */
    sort_along_vec_t *data;        /* Buffer variable needed for position sorting    */
    real             *mpi_inbuf;   /* MPI buffer                                     */
    real             *mpi_outbuf;  /* MPI buffer                                     */
    int               mpi_bufsize; /* Allocation size of in & outbuf                 */
    unsigned long     Flags;       /* mdrun flags                                    */
    gmx_bool          bOut;        /* Used to skip first output when appending to
                                    * avoid duplicate entries in rotation outfiles   */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
} t_gmx_enfrot;


/* Global enforced rotation data for a single rotation group                  */
typedef struct gmx_enfrotgrp
{
<<<<<<< HEAD
    real    degangle;       /* Rotation angle in degrees                      */
    matrix  rotmat;         /* Rotation matrix                                */
    atom_id *ind_loc;       /* Local rotation indices                         */
    int     nat_loc;        /* Number of local group atoms                    */
    int     nalloc_loc;     /* Allocation size for ind_loc and weight_loc     */

    real  V;                /* Rotation potential for this rotation group     */
    rvec  *f_rot_loc;       /* Array to store the forces on the local atoms
                               resulting from enforced rotation potential     */

    /* Collective coordinates for the whole rotation group */
    real  *xc_ref_length;   /* Length of each x_rotref vector after x_rotref 
                               has been put into origin                       */
    int   *xc_ref_ind;      /* Position of each local atom in the collective
                               array                                          */
    rvec  xc_center;        /* Center of the rotation group positions, may
                               be mass weighted                               */
    rvec  xc_ref_center;    /* dito, for the reference positions              */
=======
    real     degangle;      /* Rotation angle in degrees                      */
    matrix   rotmat;        /* Rotation matrix                                */
    atom_id *ind_loc;       /* Local rotation indices                         */
    int      nat_loc;       /* Number of local group atoms                    */
    int      nalloc_loc;    /* Allocation size for ind_loc and weight_loc     */

    real     V;             /* Rotation potential for this rotation group     */
    rvec    *f_rot_loc;     /* Array to store the forces on the local atoms
                               resulting from enforced rotation potential     */

    /* Collective coordinates for the whole rotation group */
    real  *xc_ref_length;   /* Length of each x_rotref vector after x_rotref
                               has been put into origin                       */
    int   *xc_ref_ind;      /* Position of each local atom in the collective
                               array                                          */
    rvec   xc_center;       /* Center of the rotation group positions, may
                               be mass weighted                               */
    rvec   xc_ref_center;   /* dito, for the reference positions              */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rvec  *xc;              /* Current (collective) positions                 */
    ivec  *xc_shifts;       /* Current (collective) shifts                    */
    ivec  *xc_eshifts;      /* Extra shifts since last DD step                */
    rvec  *xc_old;          /* Old (collective) positions                     */
    rvec  *xc_norm;         /* Normalized form of the current positions       */
<<<<<<< HEAD
    rvec  *xc_ref_sorted;   /* Reference positions (sorted in the same order 
=======
    rvec  *xc_ref_sorted;   /* Reference positions (sorted in the same order
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                               as xc when sorted)                             */
    int   *xc_sortind;      /* Where is a position found after sorting?       */
    real  *mc;              /* Collective masses                              */
    real  *mc_sorted;
<<<<<<< HEAD
    real  invmass;          /* one over the total mass of the rotation group  */

    real  torque_v;         /* Torque in the direction of rotation vector     */
    real  angle_v;          /* Actual angle of the whole rotation group       */
    /* Fixed rotation only */
    real  weight_v;         /* Weights for angle determination                */
=======
    real   invmass;         /* one over the total mass of the rotation group  */

    real   torque_v;        /* Torque in the direction of rotation vector     */
    real   angle_v;         /* Actual angle of the whole rotation group       */
    /* Fixed rotation only */
    real   weight_v;        /* Weights for angle determination                */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rvec  *xr_loc;          /* Local reference coords, correctly rotated      */
    rvec  *x_loc_pbc;       /* Local current coords, correct PBC image        */
    real  *m_loc;           /* Masses of the current local atoms              */

    /* Flexible rotation only */
<<<<<<< HEAD
    int   nslabs_alloc;     /* For this many slabs memory is allocated        */
    int   slab_first;       /* Lowermost slab for that the calculation needs 
                               to be performed at a given time step           */
    int   slab_last;        /* Uppermost slab ...                             */
    int   slab_first_ref;   /* First slab for which ref. center is stored     */
    int   slab_last_ref;    /* Last ...                                       */
    int   slab_buffer;      /* Slab buffer region around reference slabs      */
    int   *firstatom;       /* First relevant atom for a slab                 */
    int   *lastatom;        /* Last relevant atom for a slab                  */
    rvec  *slab_center;     /* Gaussian-weighted slab center                  */
    rvec  *slab_center_ref; /* Gaussian-weighted slab center for the
                               reference positions                            */
    real  *slab_weights;    /* Sum of gaussian weights in a slab              */
    real  *slab_torque_v;   /* Torque T = r x f for each slab.                */
                            /* torque_v = m.v = angular momentum in the 
                               direction of v                                 */
    real  max_beta;         /* min_gaussian from inputrec->rotgrp is the
                               minimum value the gaussian must have so that 
                               the force is actually evaluated max_beta is 
                               just another way to put it                     */
    real  *gn_atom;         /* Precalculated gaussians for a single atom      */
    int   *gn_slabind;      /* Tells to which slab each precalculated gaussian 
                               belongs                                        */
    rvec  *slab_innersumvec;/* Inner sum of the flexible2 potential per slab;
                               this is precalculated for optimization reasons */
    t_gmx_slabdata *slab_data; /* Holds atom positions and gaussian weights 
                               of atoms belonging to a slab                   */
=======
    int    nslabs_alloc;              /* For this many slabs memory is allocated        */
    int    slab_first;                /* Lowermost slab for that the calculation needs
                                         to be performed at a given time step           */
    int    slab_last;                 /* Uppermost slab ...                             */
    int    slab_first_ref;            /* First slab for which ref. center is stored     */
    int    slab_last_ref;             /* Last ...                                       */
    int    slab_buffer;               /* Slab buffer region around reference slabs      */
    int   *firstatom;                 /* First relevant atom for a slab                 */
    int   *lastatom;                  /* Last relevant atom for a slab                  */
    rvec  *slab_center;               /* Gaussian-weighted slab center                  */
    rvec  *slab_center_ref;           /* Gaussian-weighted slab center for the
                                         reference positions                            */
    real  *slab_weights;              /* Sum of gaussian weights in a slab              */
    real  *slab_torque_v;             /* Torque T = r x f for each slab.                */
                                      /* torque_v = m.v = angular momentum in the
                                         direction of v                                 */
    real  max_beta;                   /* min_gaussian from inputrec->rotgrp is the
                                         minimum value the gaussian must have so that
                                         the force is actually evaluated max_beta is
                                         just another way to put it                     */
    real           *gn_atom;          /* Precalculated gaussians for a single atom      */
    int            *gn_slabind;       /* Tells to which slab each precalculated gaussian
                                         belongs                                        */
    rvec           *slab_innersumvec; /* Inner sum of the flexible2 potential per slab;
                                         this is precalculated for optimization reasons */
    t_gmx_slabdata *slab_data;        /* Holds atom positions and gaussian weights
                                         of atoms belonging to a slab                   */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* For potential fits with varying angle: */
    t_gmx_potfit *PotAngleFit;  /* Used for fit type 'potential'              */
} t_gmx_enfrotgrp;


/* Activate output of forces for correctness checks */
/* #define PRINT_FORCES */
#ifdef PRINT_FORCES
<<<<<<< HEAD
#define PRINT_FORCE_J  fprintf(stderr,"f%d = %15.8f %15.8f %15.8f\n",erg->xc_ref_ind[j],erg->f_rot_loc[j][XX], erg->f_rot_loc[j][YY], erg->f_rot_loc[j][ZZ]);
#define PRINT_POT_TAU  if (MASTER(cr)) { \
                           fprintf(stderr,"potential = %15.8f\n" "torque    = %15.8f\n", erg->V, erg->torque_v); \
                       }
=======
#define PRINT_FORCE_J  fprintf(stderr, "f%d = %15.8f %15.8f %15.8f\n", erg->xc_ref_ind[j], erg->f_rot_loc[j][XX], erg->f_rot_loc[j][YY], erg->f_rot_loc[j][ZZ]);
#define PRINT_POT_TAU  if (MASTER(cr)) { \
        fprintf(stderr, "potential = %15.8f\n" "torque    = %15.8f\n", erg->V, erg->torque_v); \
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#else
#define PRINT_FORCE_J
#define PRINT_POT_TAU
#endif

/* Shortcuts for often used queries */
<<<<<<< HEAD
#define ISFLEX(rg) ( (rg->eType==erotgFLEX) || (rg->eType==erotgFLEXT) || (rg->eType==erotgFLEX2) || (rg->eType==erotgFLEX2T) )
#define ISCOLL(rg) ( (rg->eType==erotgFLEX) || (rg->eType==erotgFLEXT) || (rg->eType==erotgFLEX2) || (rg->eType==erotgFLEX2T) || (rg->eType==erotgRMPF) || (rg->eType==erotgRM2PF) )
=======
#define ISFLEX(rg) ( (rg->eType == erotgFLEX) || (rg->eType == erotgFLEXT) || (rg->eType == erotgFLEX2) || (rg->eType == erotgFLEX2T) )
#define ISCOLL(rg) ( (rg->eType == erotgFLEX) || (rg->eType == erotgFLEXT) || (rg->eType == erotgFLEX2) || (rg->eType == erotgFLEX2T) || (rg->eType == erotgRMPF) || (rg->eType == erotgRM2PF) )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/* Does any of the rotation groups use slab decomposition? */
static gmx_bool HaveFlexibleGroups(t_rot *rot)
{
<<<<<<< HEAD
    int g;
    t_rotgrp *rotg;


    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (ISFLEX(rotg))
            return TRUE;
=======
    int       g;
    t_rotgrp *rotg;


    for (g = 0; g < rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (ISFLEX(rotg))
        {
            return TRUE;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return FALSE;
}


/* Is for any group the fit angle determined by finding the minimum of the
 * rotation potential? */
static gmx_bool HavePotFitGroups(t_rot *rot)
{
<<<<<<< HEAD
    int g;
    t_rotgrp *rotg;


    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (erotgFitPOT == rotg->eFittype)
            return TRUE;
=======
    int       g;
    t_rotgrp *rotg;


    for (g = 0; g < rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (erotgFitPOT == rotg->eFittype)
        {
            return TRUE;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return FALSE;
}


static double** allocate_square_matrix(int dim)
{
<<<<<<< HEAD
    int i;
    double** mat = NULL; 
    
    
    snew(mat, dim);
    for(i=0; i<dim; i++)
        snew(mat[i], dim);
=======
    int      i;
    double** mat = NULL;


    snew(mat, dim);
    for (i = 0; i < dim; i++)
    {
        snew(mat[i], dim);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return mat;
}


static void free_square_matrix(double** mat, int dim)
{
    int i;
<<<<<<< HEAD
    
    
    for (i=0; i<dim; i++)
        sfree(mat[i]);
=======


    for (i = 0; i < dim; i++)
    {
        sfree(mat[i]);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    sfree(mat);
}


/* Return the angle for which the potential is minimal */
static real get_fitangle(t_rotgrp *rotg, gmx_enfrotgrp_t erg)
{
<<<<<<< HEAD
    int i;
    real fitangle = -999.9;
    real pot_min = GMX_FLOAT_MAX;
=======
    int           i;
    real          fitangle = -999.9;
    real          pot_min  = GMX_FLOAT_MAX;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_gmx_potfit *fit;


    fit = erg->PotAngleFit;

    for (i = 0; i < rotg->PotAngle_nstep; i++)
    {
        if (fit->V[i] < pot_min)
        {
<<<<<<< HEAD
            pot_min = fit->V[i];
=======
            pot_min  = fit->V[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            fitangle = fit->degangle[i];
        }
    }

    return fitangle;
}


/* Reduce potential angle fit data for this group at this time step? */
static gmx_inline gmx_bool bPotAngle(t_rot *rot, t_rotgrp *rotg, gmx_large_int_t step)
{
<<<<<<< HEAD
    return ( (erotgFitPOT==rotg->eFittype) && (do_per_step(step, rot->nstsout) || do_per_step(step, rot->nstrout)) );
=======
    return ( (erotgFitPOT == rotg->eFittype) && (do_per_step(step, rot->nstsout) || do_per_step(step, rot->nstrout)) );
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

/* Reduce slab torqe data for this group at this time step? */
static gmx_inline gmx_bool bSlabTau(t_rot *rot, t_rotgrp *rotg, gmx_large_int_t step)
{
    return ( (ISFLEX(rotg)) && do_per_step(step, rot->nstsout) );
}

/* Output rotation energy, torques, etc. for each rotation group */
static void reduce_output(t_commrec *cr, t_rot *rot, real t, gmx_large_int_t step)
{
<<<<<<< HEAD
    int      g,i,islab,nslabs=0;
    int      count;      /* MPI element counter                               */
    t_rotgrp *rotg;
    gmx_enfrot_t er;     /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg; /* Pointer to enforced rotation group data           */
    real     fitangle;
    gmx_bool bFlex;

    
    er=rot->enfrot;
    
=======
    int             g, i, islab, nslabs = 0;
    int             count; /* MPI element counter                               */
    t_rotgrp       *rotg;
    gmx_enfrot_t    er;    /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg;   /* Pointer to enforced rotation group data           */
    real            fitangle;
    gmx_bool        bFlex;


    er = rot->enfrot;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Fill the MPI buffer with stuff to reduce. If items are added for reduction
     * here, the MPI buffer size has to be enlarged also in calc_mpi_bufsize() */
    if (PAR(cr))
    {
<<<<<<< HEAD
        count=0;
        for (g=0; g < rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            erg = rotg->enfrotgrp;
            nslabs = erg->slab_last - erg->slab_first + 1;
=======
        count = 0;
        for (g = 0; g < rot->ngrp; g++)
        {
            rotg                   = &rot->grp[g];
            erg                    = rotg->enfrotgrp;
            nslabs                 = erg->slab_last - erg->slab_first + 1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            er->mpi_inbuf[count++] = erg->V;
            er->mpi_inbuf[count++] = erg->torque_v;
            er->mpi_inbuf[count++] = erg->angle_v;
            er->mpi_inbuf[count++] = erg->weight_v; /* weights are not needed for flex types, but this is just a single value */

            if (bPotAngle(rot, rotg, step))
            {
                for (i = 0; i < rotg->PotAngle_nstep; i++)
<<<<<<< HEAD
                    er->mpi_inbuf[count++] = erg->PotAngleFit->V[i];
            }
            if (bSlabTau(rot, rotg, step))
            {
                for (i=0; i<nslabs; i++)
                    er->mpi_inbuf[count++] = erg->slab_torque_v[i];
            }
        }
        if (count > er->mpi_bufsize)
            gmx_fatal(FARGS, "%s MPI buffer overflow, please report this error.", RotStr);
=======
                {
                    er->mpi_inbuf[count++] = erg->PotAngleFit->V[i];
                }
            }
            if (bSlabTau(rot, rotg, step))
            {
                for (i = 0; i < nslabs; i++)
                {
                    er->mpi_inbuf[count++] = erg->slab_torque_v[i];
                }
            }
        }
        if (count > er->mpi_bufsize)
        {
            gmx_fatal(FARGS, "%s MPI buffer overflow, please report this error.", RotStr);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef GMX_MPI
        MPI_Reduce(er->mpi_inbuf, er->mpi_outbuf, count, GMX_MPI_REAL, MPI_SUM, MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif

        /* Copy back the reduced data from the buffer on the master */
        if (MASTER(cr))
        {
<<<<<<< HEAD
            count=0;
            for (g=0; g < rot->ngrp; g++)
            {
                rotg = &rot->grp[g];
                erg = rotg->enfrotgrp;
                nslabs = erg->slab_last - erg->slab_first + 1;
=======
            count = 0;
            for (g = 0; g < rot->ngrp; g++)
            {
                rotg          = &rot->grp[g];
                erg           = rotg->enfrotgrp;
                nslabs        = erg->slab_last - erg->slab_first + 1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                erg->V        = er->mpi_outbuf[count++];
                erg->torque_v = er->mpi_outbuf[count++];
                erg->angle_v  = er->mpi_outbuf[count++];
                erg->weight_v = er->mpi_outbuf[count++];

                if (bPotAngle(rot, rotg, step))
                {
                    for (i = 0; i < rotg->PotAngle_nstep; i++)
<<<<<<< HEAD
                        erg->PotAngleFit->V[i] = er->mpi_outbuf[count++];
                }
                if (bSlabTau(rot, rotg, step))
                {
                    for (i=0; i<nslabs; i++)
                        erg->slab_torque_v[i] = er->mpi_outbuf[count++];
=======
                    {
                        erg->PotAngleFit->V[i] = er->mpi_outbuf[count++];
                    }
                }
                if (bSlabTau(rot, rotg, step))
                {
                    for (i = 0; i < nslabs; i++)
                    {
                        erg->slab_torque_v[i] = er->mpi_outbuf[count++];
                    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
            }
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Output */
    if (MASTER(cr))
    {
        /* Angle and torque for each rotation group */
<<<<<<< HEAD
        for (g=0; g < rot->ngrp; g++)
        {
            rotg=&rot->grp[g];
            bFlex = ISFLEX(rotg);

            erg=rotg->enfrotgrp;
            
            /* Output to main rotation output file: */
            if ( do_per_step(step, rot->nstrout) )
=======
        for (g = 0; g < rot->ngrp; g++)
        {
            rotg  = &rot->grp[g];
            bFlex = ISFLEX(rotg);

            erg = rotg->enfrotgrp;

            /* Output to main rotation output file: */
            if (do_per_step(step, rot->nstrout) )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                if (erotgFitPOT == rotg->eFittype)
                {
                    fitangle = get_fitangle(rotg, erg);
                }
                else
                {
                    if (bFlex)
<<<<<<< HEAD
                        fitangle = erg->angle_v; /* RMSD fit angle */
                    else
                        fitangle = (erg->angle_v/erg->weight_v)*180.0*M_1_PI;
=======
                    {
                        fitangle = erg->angle_v; /* RMSD fit angle */
                    }
                    else
                    {
                        fitangle = (erg->angle_v/erg->weight_v)*180.0*M_1_PI;
                    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
                fprintf(er->out_rot, "%12.4f", fitangle);
                fprintf(er->out_rot, "%12.3e", erg->torque_v);
                fprintf(er->out_rot, "%12.3e", erg->V);
            }

<<<<<<< HEAD
            if ( do_per_step(step, rot->nstsout) )
=======
            if (do_per_step(step, rot->nstsout) )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                /* Output to torque log file: */
                if (bFlex)
                {
                    fprintf(er->out_torque, "%12.3e%6d", t, g);
<<<<<<< HEAD
                    for (i=erg->slab_first; i<=erg->slab_last; i++)
=======
                    for (i = erg->slab_first; i <= erg->slab_last; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        islab = i - erg->slab_first;  /* slab index */
                        /* Only output if enough weight is in slab */
                        if (erg->slab_weights[islab] > rotg->min_gaussian)
<<<<<<< HEAD
                            fprintf(er->out_torque, "%6d%12.3e", i, erg->slab_torque_v[islab]);
                    }
                    fprintf(er->out_torque , "\n");
=======
                        {
                            fprintf(er->out_torque, "%6d%12.3e", i, erg->slab_torque_v[islab]);
                        }
                    }
                    fprintf(er->out_torque, "\n");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }

                /* Output to angles log file: */
                if (erotgFitPOT == rotg->eFittype)
                {
                    fprintf(er->out_angles, "%12.3e%6d%12.4f", t, g, erg->degangle);
                    /* Output energies at a set of angles around the reference angle */
                    for (i = 0; i < rotg->PotAngle_nstep; i++)
<<<<<<< HEAD
                        fprintf(er->out_angles, "%12.3e", erg->PotAngleFit->V[i]);
=======
                    {
                        fprintf(er->out_angles, "%12.3e", erg->PotAngleFit->V[i]);
                    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    fprintf(er->out_angles, "\n");
                }
            }
        }
<<<<<<< HEAD
        if ( do_per_step(step, rot->nstrout) )
            fprintf(er->out_rot, "\n");
=======
        if (do_per_step(step, rot->nstrout) )
        {
            fprintf(er->out_rot, "\n");
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}


/* Add the forces from enforced rotation potential to the local forces.
 * Should be called after the SR forces have been evaluated */
extern real add_rot_forces(t_rot *rot, rvec f[], t_commrec *cr, gmx_large_int_t step, real t)
{
<<<<<<< HEAD
    int g,l,ii;
    t_rotgrp *rotg;
    gmx_enfrot_t er;     /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg; /* Pointer to enforced rotation group data           */
    real Vrot = 0.0;     /* If more than one rotation group is present, Vrot
                            assembles the local parts from all groups         */

    
    er=rot->enfrot;
    
    /* Loop over enforced rotation groups (usually 1, though)
     * Apply the forces from rotation potentials */
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg=rotg->enfrotgrp;
        Vrot += erg->V;  /* add the local parts from the nodes */
        for (l=0; l<erg->nat_loc; l++)
=======
    int             g, l, ii;
    t_rotgrp       *rotg;
    gmx_enfrot_t    er;         /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data           */
    real            Vrot = 0.0; /* If more than one rotation group is present, Vrot
                                   assembles the local parts from all groups         */


    er = rot->enfrot;

    /* Loop over enforced rotation groups (usually 1, though)
     * Apply the forces from rotation potentials */
    for (g = 0; g < rot->ngrp; g++)
    {
        rotg  = &rot->grp[g];
        erg   = rotg->enfrotgrp;
        Vrot += erg->V;  /* add the local parts from the nodes */
        for (l = 0; l < erg->nat_loc; l++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* Get the right index of the local force */
            ii = erg->ind_loc[l];
            /* Add */
<<<<<<< HEAD
            rvec_inc(f[ii],erg->f_rot_loc[l]);
=======
            rvec_inc(f[ii], erg->f_rot_loc[l]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }

    /* Reduce energy,torque, angles etc. to get the sum values (per rotation group)
     * on the master and output these values to file. */
    if ( (do_per_step(step, rot->nstrout) || do_per_step(step, rot->nstsout)) && er->bOut)
<<<<<<< HEAD
        reduce_output(cr, rot, t, step);

    /* When appending, er->bOut is FALSE the first time to avoid duplicate entries */
    er->bOut = TRUE;
    
=======
    {
        reduce_output(cr, rot, t, step);
    }

    /* When appending, er->bOut is FALSE the first time to avoid duplicate entries */
    er->bOut = TRUE;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    PRINT_POT_TAU

    return Vrot;
}


/* The Gaussian norm is chosen such that the sum of the gaussian functions
 * over the slabs is approximately 1.0 everywhere */
#define GAUSS_NORM   0.569917543430618


/* Calculate the maximum beta that leads to a gaussian larger min_gaussian,
 * also does some checks
 */
static double calc_beta_max(real min_gaussian, real slab_dist)
{
    double sigma;
    double arg;
<<<<<<< HEAD
    
    
    /* Actually the next two checks are already made in grompp */
    if (slab_dist <= 0)
        gmx_fatal(FARGS, "Slab distance of flexible rotation groups must be >=0 !");
    if (min_gaussian <= 0)
        gmx_fatal(FARGS, "Cutoff value for Gaussian must be > 0. (You requested %f)");
=======


    /* Actually the next two checks are already made in grompp */
    if (slab_dist <= 0)
    {
        gmx_fatal(FARGS, "Slab distance of flexible rotation groups must be >=0 !");
    }
    if (min_gaussian <= 0)
    {
        gmx_fatal(FARGS, "Cutoff value for Gaussian must be > 0. (You requested %f)");
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Define the sigma value */
    sigma = 0.7*slab_dist;

    /* Calculate the argument for the logarithm and check that the log() result is negative or 0 */
    arg = min_gaussian/GAUSS_NORM;
    if (arg > 1.0)
<<<<<<< HEAD
        gmx_fatal(FARGS, "min_gaussian of flexible rotation groups must be <%g", GAUSS_NORM);
    
=======
    {
        gmx_fatal(FARGS, "min_gaussian of flexible rotation groups must be <%g", GAUSS_NORM);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return sqrt(-2.0*sigma*sigma*log(min_gaussian/GAUSS_NORM));
}


static gmx_inline real calc_beta(rvec curr_x, t_rotgrp *rotg, int n)
{
    return iprod(curr_x, rotg->vec) - rotg->slab_dist * n;
}


static gmx_inline real gaussian_weight(rvec curr_x, t_rotgrp *rotg, int n)
{
    const real norm = GAUSS_NORM;
    real       sigma;

<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;
    /* Calculate the Gaussian value of slab n for position curr_x */
    return norm * exp( -0.5 * sqr( calc_beta(curr_x, rotg, n)/sigma ) );
}


/* Returns the weight in a single slab, also calculates the Gaussian- and mass-
 * weighted sum of positions for that slab */
static real get_slab_weight(int j, t_rotgrp *rotg, rvec xc[], real mc[], rvec *x_weighted_sum)
{
<<<<<<< HEAD
    rvec curr_x;              /* The position of an atom                      */
    rvec curr_x_weighted;     /* The gaussian-weighted position               */
    real gaussian;            /* A single gaussian weight                     */
    real wgauss;              /* gaussian times current mass                  */
    real slabweight = 0.0;    /* The sum of weights in the slab               */
    int i,islab;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data      */

    
    erg=rotg->enfrotgrp;
    clear_rvec(*x_weighted_sum);
    
    /* Slab index */
    islab = j - erg->slab_first;
    
    /* Loop over all atoms in the rotation group */
     for (i=0; i<rotg->nat; i++)
     {
         copy_rvec(xc[i], curr_x);
         gaussian = gaussian_weight(curr_x, rotg, j);
         wgauss = gaussian * mc[i];
         svmul(wgauss, curr_x, curr_x_weighted);
         rvec_add(*x_weighted_sum, curr_x_weighted, *x_weighted_sum);
         slabweight += wgauss;
     } /* END of loop over rotation group atoms */

     return slabweight;
=======
    rvec            curr_x;           /* The position of an atom                      */
    rvec            curr_x_weighted;  /* The gaussian-weighted position               */
    real            gaussian;         /* A single gaussian weight                     */
    real            wgauss;           /* gaussian times current mass                  */
    real            slabweight = 0.0; /* The sum of weights in the slab               */
    int             i, islab;
    gmx_enfrotgrp_t erg;              /* Pointer to enforced rotation group data      */


    erg = rotg->enfrotgrp;
    clear_rvec(*x_weighted_sum);

    /* Slab index */
    islab = j - erg->slab_first;

    /* Loop over all atoms in the rotation group */
    for (i = 0; i < rotg->nat; i++)
    {
        copy_rvec(xc[i], curr_x);
        gaussian = gaussian_weight(curr_x, rotg, j);
        wgauss   = gaussian * mc[i];
        svmul(wgauss, curr_x, curr_x_weighted);
        rvec_add(*x_weighted_sum, curr_x_weighted, *x_weighted_sum);
        slabweight += wgauss;
    }  /* END of loop over rotation group atoms */

    return slabweight;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static void get_slab_centers(
<<<<<<< HEAD
        t_rotgrp *rotg,       /* The rotation group information               */
        rvec      *xc,        /* The rotation group positions; will 
                                 typically be enfrotgrp->xc, but at first call 
                                 it is enfrotgrp->xc_ref                      */
        real      *mc,        /* The masses of the rotation group atoms       */
        int       g,          /* The number of the rotation group             */
        real      time,       /* Used for output only                         */
        FILE      *out_slabs, /* For outputting center per slab information   */
        gmx_bool  bOutStep,   /* Is this an output step?                      */
        gmx_bool  bReference) /* If this routine is called from
                                 init_rot_group we need to store
                                 the reference slab centers                   */
{
    int j,islab;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    
    
    erg=rotg->enfrotgrp;
=======
        t_rotgrp  *rotg,       /* The rotation group information               */
        rvec      *xc,         /* The rotation group positions; will
                                  typically be enfrotgrp->xc, but at first call
                                  it is enfrotgrp->xc_ref                      */
        real      *mc,         /* The masses of the rotation group atoms       */
        int        g,          /* The number of the rotation group             */
        real       time,       /* Used for output only                         */
        FILE      *out_slabs,  /* For outputting center per slab information   */
        gmx_bool   bOutStep,   /* Is this an output step?                      */
        gmx_bool   bReference) /* If this routine is called from
                                  init_rot_group we need to store
                                  the reference slab centers                   */
{
    int             j, islab;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Loop over slabs */
    for (j = erg->slab_first; j <= erg->slab_last; j++)
    {
<<<<<<< HEAD
        islab = j - erg->slab_first;
        erg->slab_weights[islab] = get_slab_weight(j, rotg, xc, mc, &erg->slab_center[islab]);
        
=======
        islab                    = j - erg->slab_first;
        erg->slab_weights[islab] = get_slab_weight(j, rotg, xc, mc, &erg->slab_center[islab]);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* We can do the calculations ONLY if there is weight in the slab! */
        if (erg->slab_weights[islab] > WEIGHT_MIN)
        {
            svmul(1.0/erg->slab_weights[islab], erg->slab_center[islab], erg->slab_center[islab]);
        }
        else
        {
            /* We need to check this here, since we divide through slab_weights
             * in the flexible low-level routines! */
            gmx_fatal(FARGS, "Not enough weight in slab %d. Slab center cannot be determined!", j);
        }
<<<<<<< HEAD
        
        /* At first time step: save the centers of the reference structure */
        if (bReference)
            copy_rvec(erg->slab_center[islab], erg->slab_center_ref[islab]);
    } /* END of loop over slabs */
    
=======

        /* At first time step: save the centers of the reference structure */
        if (bReference)
        {
            copy_rvec(erg->slab_center[islab], erg->slab_center_ref[islab]);
        }
    } /* END of loop over slabs */

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Output on the master */
    if ( (NULL != out_slabs) && bOutStep)
    {
        fprintf(out_slabs, "%12.3e%6d", time, g);
        for (j = erg->slab_first; j <= erg->slab_last; j++)
        {
            islab = j - erg->slab_first;
            fprintf(out_slabs, "%6d%12.3e%12.3e%12.3e",
<<<<<<< HEAD
                    j,erg->slab_center[islab][XX],erg->slab_center[islab][YY],erg->slab_center[islab][ZZ]);
=======
                    j, erg->slab_center[islab][XX], erg->slab_center[islab][YY], erg->slab_center[islab][ZZ]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        fprintf(out_slabs, "\n");
    }
}


static void calc_rotmat(
<<<<<<< HEAD
        rvec vec,
        real degangle,  /* Angle alpha of rotation at time t in degrees       */
        matrix rotmat)  /* Rotation matrix                                    */
=======
        rvec   vec,
        real   degangle,      /* Angle alpha of rotation at time t in degrees       */
        matrix rotmat)        /* Rotation matrix                                    */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    real radangle;            /* Rotation angle in radians */
    real cosa;                /* cosine alpha              */
    real sina;                /* sine alpha                */
    real OMcosa;              /* 1 - cos(alpha)            */
    real dumxy, dumxz, dumyz; /* save computations         */
    rvec rot_vec;             /* Rotate around rot_vec ... */


    radangle = degangle * M_PI/180.0;
<<<<<<< HEAD
    copy_rvec(vec , rot_vec );
=======
    copy_rvec(vec, rot_vec );
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Precompute some variables: */
    cosa   = cos(radangle);
    sina   = sin(radangle);
    OMcosa = 1.0 - cosa;
    dumxy  = rot_vec[XX]*rot_vec[YY]*OMcosa;
    dumxz  = rot_vec[XX]*rot_vec[ZZ]*OMcosa;
    dumyz  = rot_vec[YY]*rot_vec[ZZ]*OMcosa;

    /* Construct the rotation matrix for this rotation group: */
    /* 1st column: */
    rotmat[XX][XX] = cosa  + rot_vec[XX]*rot_vec[XX]*OMcosa;
    rotmat[YY][XX] = dumxy + rot_vec[ZZ]*sina;
    rotmat[ZZ][XX] = dumxz - rot_vec[YY]*sina;
    /* 2nd column: */
    rotmat[XX][YY] = dumxy - rot_vec[ZZ]*sina;
    rotmat[YY][YY] = cosa  + rot_vec[YY]*rot_vec[YY]*OMcosa;
    rotmat[ZZ][YY] = dumyz + rot_vec[XX]*sina;
    /* 3rd column: */
    rotmat[XX][ZZ] = dumxz + rot_vec[YY]*sina;
    rotmat[YY][ZZ] = dumyz - rot_vec[XX]*sina;
    rotmat[ZZ][ZZ] = cosa  + rot_vec[ZZ]*rot_vec[ZZ]*OMcosa;

#ifdef PRINTMATRIX
<<<<<<< HEAD
    int iii,jjj;

    for (iii=0; iii<3; iii++) {
        for (jjj=0; jjj<3; jjj++)
            fprintf(stderr, " %10.8f ",  rotmat[iii][jjj]);
=======
    int iii, jjj;

    for (iii = 0; iii < 3; iii++)
    {
        for (jjj = 0; jjj < 3; jjj++)
        {
            fprintf(stderr, " %10.8f ",  rotmat[iii][jjj]);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fprintf(stderr, "\n");
    }
#endif
}


/* Calculates torque on the rotation axis tau = position x force */
static gmx_inline real torque(
        rvec rotvec,  /* rotation vector; MUST be normalized!                 */
        rvec force,   /* force                                                */
        rvec x,       /* position of atom on which the force acts             */
        rvec pivot)   /* pivot point of rotation axis                         */
{
    rvec vectmp, tau;

<<<<<<< HEAD
    
    /* Subtract offset */
    rvec_sub(x,pivot,vectmp);
    
    /* position x force */
    cprod(vectmp, force, tau);
    
=======

    /* Subtract offset */
    rvec_sub(x, pivot, vectmp);

    /* position x force */
    cprod(vectmp, force, tau);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Return the part of the torque which is parallel to the rotation vector */
    return iprod(tau, rotvec);
}


/* Right-aligned output of value with standard width */
static void print_aligned(FILE *fp, char *str)
{
    fprintf(fp, "%12s", str);
}


/* Right-aligned output of value with standard short width */
static void print_aligned_short(FILE *fp, char *str)
{
    fprintf(fp, "%6s", str);
}


static FILE *open_output_file(const char *fn, int steps, const char what[])
{
    FILE *fp;
<<<<<<< HEAD
    
    
    fp = ffopen(fn, "w");

    fprintf(fp, "# Output of %s is written in intervals of %d time step%s.\n#\n",
            what,steps, steps>1 ? "s":"");
    
=======


    fp = ffopen(fn, "w");

    fprintf(fp, "# Output of %s is written in intervals of %d time step%s.\n#\n",
            what, steps, steps > 1 ? "s" : "");

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return fp;
}


/* Open output file for slab center data. Call on master only */
static FILE *open_slab_out(const char *fn, t_rot *rot, const output_env_t oenv)
{
    FILE      *fp;
<<<<<<< HEAD
    int       g,i;
=======
    int        g, i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_rotgrp  *rotg;


    if (rot->enfrot->Flags & MD_APPENDFILES)
    {
<<<<<<< HEAD
        fp = gmx_fio_fopen(fn,"a");
=======
        fp = gmx_fio_fopen(fn, "a");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        fp = open_output_file(fn, rot->nstsout, "gaussian weighted slab centers");

<<<<<<< HEAD
        for (g=0; g<rot->ngrp; g++)
=======
        for (g = 0; g < rot->ngrp; g++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            rotg = &rot->grp[g];
            if (ISFLEX(rotg))
            {
                fprintf(fp, "# Rotation group %d (%s), slab distance %f nm, %s.\n",
                        g, erotg_names[rotg->eType], rotg->slab_dist,
<<<<<<< HEAD
                        rotg->bMassW? "centers of mass":"geometrical centers");
=======
                        rotg->bMassW ? "centers of mass" : "geometrical centers");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }

        fprintf(fp, "# Reference centers are listed first (t=-1).\n");
        fprintf(fp, "# The following columns have the syntax:\n");
        fprintf(fp, "#     ");
        print_aligned_short(fp, "t");
        print_aligned_short(fp, "grp");
        /* Print legend for the first two entries only ... */
<<<<<<< HEAD
        for (i=0; i<2; i++)
=======
        for (i = 0; i < 2; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            print_aligned_short(fp, "slab");
            print_aligned(fp, "X center");
            print_aligned(fp, "Y center");
            print_aligned(fp, "Z center");
        }
        fprintf(fp, " ...\n");
        fflush(fp);
    }

    return fp;
}


/* Adds 'buf' to 'str' */
static void add_to_string(char **str, char *buf)
{
    int len;


    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}


static void add_to_string_aligned(char **str, char *buf)
{
    char buf_aligned[STRLEN];

    sprintf(buf_aligned, "%12s", buf);
    add_to_string(str, buf_aligned);
}


/* Open output file and print some general information about the rotation groups.
 * Call on master only */
static FILE *open_rot_out(const char *fn, t_rot *rot, const output_env_t oenv)
{
<<<<<<< HEAD
    FILE       *fp;
    int        g,nsets;
    t_rotgrp   *rotg;
    const char **setname;
    char       buf[50], buf2[75];
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    gmx_bool   bFlex;
    char       *LegendStr=NULL;
=======
    FILE           *fp;
    int             g, nsets;
    t_rotgrp       *rotg;
    const char    **setname;
    char            buf[50], buf2[75];
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    gmx_bool        bFlex;
    char           *LegendStr = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    if (rot->enfrot->Flags & MD_APPENDFILES)
    {
<<<<<<< HEAD
        fp = gmx_fio_fopen(fn,"a");
=======
        fp = gmx_fio_fopen(fn, "a");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        fp = xvgropen(fn, "Rotation angles and energy", "Time (ps)", "angles (degrees) and energies (kJ/mol)", oenv);
<<<<<<< HEAD
        fprintf(fp, "# Output of enforced rotation data is written in intervals of %d time step%s.\n#\n", rot->nstrout, rot->nstrout > 1 ? "s":"");
=======
        fprintf(fp, "# Output of enforced rotation data is written in intervals of %d time step%s.\n#\n", rot->nstrout, rot->nstrout > 1 ? "s" : "");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fprintf(fp, "# The scalar tau is the torque (kJ/mol) in the direction of the rotation vector v.\n");
        fprintf(fp, "# To obtain the vectorial torque, multiply tau with the group's rot_vec.\n");
        fprintf(fp, "# For flexible groups, tau(t,n) from all slabs n have been summed in a single value tau(t) here.\n");
        fprintf(fp, "# The torques tau(t,n) are found in the rottorque.log (-rt) output file\n");
<<<<<<< HEAD
        
        for (g=0; g<rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            erg=rotg->enfrotgrp;
            bFlex = ISFLEX(rotg);

            fprintf(fp, "#\n");
            fprintf(fp, "# ROTATION GROUP %d, potential type '%s':\n"      , g, erotg_names[rotg->eType]);
            fprintf(fp, "# rot_massw%d          %s\n"                      , g, yesno_names[rotg->bMassW]);
            fprintf(fp, "# rot_vec%d            %12.5e %12.5e %12.5e\n"    , g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
            fprintf(fp, "# rot_rate%d           %12.5e degrees/ps\n"       , g, rotg->rate);
            fprintf(fp, "# rot_k%d              %12.5e kJ/(mol*nm^2)\n"    , g, rotg->k);
            if ( rotg->eType==erotgISO || rotg->eType==erotgPM || rotg->eType==erotgRM || rotg->eType==erotgRM2)
                fprintf(fp, "# rot_pivot%d          %12.5e %12.5e %12.5e  nm\n", g, rotg->pivot[XX], rotg->pivot[YY], rotg->pivot[ZZ]);
=======

        for (g = 0; g < rot->ngrp; g++)
        {
            rotg  = &rot->grp[g];
            erg   = rotg->enfrotgrp;
            bFlex = ISFLEX(rotg);

            fprintf(fp, "#\n");
            fprintf(fp, "# ROTATION GROUP %d, potential type '%s':\n", g, erotg_names[rotg->eType]);
            fprintf(fp, "# rot_massw%d          %s\n", g, yesno_names[rotg->bMassW]);
            fprintf(fp, "# rot_vec%d            %12.5e %12.5e %12.5e\n", g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
            fprintf(fp, "# rot_rate%d           %12.5e degrees/ps\n", g, rotg->rate);
            fprintf(fp, "# rot_k%d              %12.5e kJ/(mol*nm^2)\n", g, rotg->k);
            if (rotg->eType == erotgISO || rotg->eType == erotgPM || rotg->eType == erotgRM || rotg->eType == erotgRM2)
            {
                fprintf(fp, "# rot_pivot%d          %12.5e %12.5e %12.5e  nm\n", g, rotg->pivot[XX], rotg->pivot[YY], rotg->pivot[ZZ]);
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            if (bFlex)
            {
                fprintf(fp, "# rot_slab_distance%d   %f nm\n", g, rotg->slab_dist);
                fprintf(fp, "# rot_min_gaussian%d   %12.5e\n", g, rotg->min_gaussian);
            }

            /* Output the centers of the rotation groups for the pivot-free potentials */
<<<<<<< HEAD
            if ((rotg->eType==erotgISOPF) || (rotg->eType==erotgPMPF) || (rotg->eType==erotgRMPF) || (rotg->eType==erotgRM2PF
                || (rotg->eType==erotgFLEXT) || (rotg->eType==erotgFLEX2T)) )
            {
                fprintf(fp, "# ref. grp. %d center  %12.5e %12.5e %12.5e\n", g,
                            erg->xc_ref_center[XX], erg->xc_ref_center[YY], erg->xc_ref_center[ZZ]);

                fprintf(fp, "# grp. %d init.center  %12.5e %12.5e %12.5e\n", g,
                            erg->xc_center[XX], erg->xc_center[YY], erg->xc_center[ZZ]);
            }

            if ( (rotg->eType == erotgRM2) || (rotg->eType==erotgFLEX2) || (rotg->eType==erotgFLEX2T) )
=======
            if ((rotg->eType == erotgISOPF) || (rotg->eType == erotgPMPF) || (rotg->eType == erotgRMPF) || (rotg->eType == erotgRM2PF
                                                                                                            || (rotg->eType == erotgFLEXT) || (rotg->eType == erotgFLEX2T)) )
            {
                fprintf(fp, "# ref. grp. %d center  %12.5e %12.5e %12.5e\n", g,
                        erg->xc_ref_center[XX], erg->xc_ref_center[YY], erg->xc_ref_center[ZZ]);

                fprintf(fp, "# grp. %d init.center  %12.5e %12.5e %12.5e\n", g,
                        erg->xc_center[XX], erg->xc_center[YY], erg->xc_center[ZZ]);
            }

            if ( (rotg->eType == erotgRM2) || (rotg->eType == erotgFLEX2) || (rotg->eType == erotgFLEX2T) )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                fprintf(fp, "# rot_eps%d            %12.5e nm^2\n", g, rotg->eps);
            }
            if (erotgFitPOT == rotg->eFittype)
            {
                fprintf(fp, "#\n");
                fprintf(fp, "# theta_fit%d is determined by first evaluating the potential for %d angles around theta_ref%d.\n",
<<<<<<< HEAD
                            g, rotg->PotAngle_nstep, g);
=======
                        g, rotg->PotAngle_nstep, g);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                fprintf(fp, "# The fit angle is the one with the smallest potential. It is given as the deviation\n");
                fprintf(fp, "# from the reference angle, i.e. if theta_ref=X and theta_fit=Y, then the angle with\n");
                fprintf(fp, "# minimal value of the potential is X+Y. Angular resolution is %g degrees.\n", rotg->PotAngle_step);
            }
        }
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Print a nice legend */
        snew(LegendStr, 1);
        LegendStr[0] = '\0';
        sprintf(buf, "#     %6s", "time");
        add_to_string_aligned(&LegendStr, buf);

        nsets = 0;
        snew(setname, 4*rot->ngrp);
<<<<<<< HEAD
        
        for (g=0; g<rot->ngrp; g++)
=======

        for (g = 0; g < rot->ngrp; g++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            rotg = &rot->grp[g];
            sprintf(buf, "theta_ref%d", g);
            add_to_string_aligned(&LegendStr, buf);

            sprintf(buf2, "%s (degrees)", buf);
            setname[nsets] = strdup(buf2);
            nsets++;
        }
<<<<<<< HEAD
        for (g=0; g<rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
=======
        for (g = 0; g < rot->ngrp; g++)
        {
            rotg  = &rot->grp[g];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            bFlex = ISFLEX(rotg);

            /* For flexible axis rotation we use RMSD fitting to determine the
             * actual angle of the rotation group */
            if (bFlex || erotgFitPOT == rotg->eFittype)
<<<<<<< HEAD
                sprintf(buf, "theta_fit%d", g);
            else
                sprintf(buf, "theta_av%d", g);
=======
            {
                sprintf(buf, "theta_fit%d", g);
            }
            else
            {
                sprintf(buf, "theta_av%d", g);
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            add_to_string_aligned(&LegendStr, buf);
            sprintf(buf2, "%s (degrees)", buf);
            setname[nsets] = strdup(buf2);
            nsets++;

            sprintf(buf, "tau%d", g);
            add_to_string_aligned(&LegendStr, buf);
            sprintf(buf2, "%s (kJ/mol)", buf);
            setname[nsets] = strdup(buf2);
            nsets++;

            sprintf(buf, "energy%d", g);
            add_to_string_aligned(&LegendStr, buf);
            sprintf(buf2, "%s (kJ/mol)", buf);
            setname[nsets] = strdup(buf2);
            nsets++;
        }
        fprintf(fp, "#\n");
<<<<<<< HEAD
        
        if (nsets > 1)
            xvgr_legend(fp, nsets, setname, oenv);
=======

        if (nsets > 1)
        {
            xvgr_legend(fp, nsets, setname, oenv);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        sfree(setname);

        fprintf(fp, "#\n# Legend for the following data columns:\n");
        fprintf(fp, "%s\n", LegendStr);
        sfree(LegendStr);
<<<<<<< HEAD
        
        fflush(fp);
    }
    
=======

        fflush(fp);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return fp;
}


/* Call on master only */
static FILE *open_angles_out(const char *fn, t_rot *rot, const output_env_t oenv)
{
<<<<<<< HEAD
    int      g,i;
    FILE     *fp;
    t_rotgrp *rotg;
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data */
    char     buf[100];
=======
    int             g, i;
    FILE           *fp;
    t_rotgrp       *rotg;
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data */
    char            buf[100];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    if (rot->enfrot->Flags & MD_APPENDFILES)
    {
<<<<<<< HEAD
        fp = gmx_fio_fopen(fn,"a");
=======
        fp = gmx_fio_fopen(fn, "a");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        /* Open output file and write some information about it's structure: */
        fp = open_output_file(fn, rot->nstsout, "rotation group angles");
        fprintf(fp, "# All angles given in degrees, time in ps.\n");
<<<<<<< HEAD
        for (g=0; g<rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            erg=rotg->enfrotgrp;

            /* Output for this group happens only if potential type is flexible or
             * if fit type is potential! */
            if ( ISFLEX(rotg) || (erotgFitPOT == rotg->eFittype) )
            {
                if (ISFLEX(rotg))
                    sprintf(buf, " slab distance %f nm, ", rotg->slab_dist);
                else
                    buf[0] = '\0';
=======
        for (g = 0; g < rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            erg  = rotg->enfrotgrp;

            /* Output for this group happens only if potential type is flexible or
             * if fit type is potential! */
            if (ISFLEX(rotg) || (erotgFitPOT == rotg->eFittype) )
            {
                if (ISFLEX(rotg))
                {
                    sprintf(buf, " slab distance %f nm, ", rotg->slab_dist);
                }
                else
                {
                    buf[0] = '\0';
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

                fprintf(fp, "#\n# ROTATION GROUP %d '%s',%s fit type '%s'.\n",
                        g, erotg_names[rotg->eType], buf, erotg_fitnames[rotg->eFittype]);

                /* Special type of fitting using the potential minimum. This is
                 * done for the whole group only, not for the individual slabs. */
                if (erotgFitPOT == rotg->eFittype)
                {
                    fprintf(fp, "#    To obtain theta_fit%d, the potential is evaluated for %d angles around theta_ref%d\n", g, rotg->PotAngle_nstep, g);
                    fprintf(fp, "#    The fit angle in the rotation standard outfile is the one with minimal energy E(theta_fit) [kJ/mol].\n");
                    fprintf(fp, "#\n");
                }

                fprintf(fp, "# Legend for the group %d data columns:\n", g);
                fprintf(fp, "#     ");
                print_aligned_short(fp, "time");
                print_aligned_short(fp, "grp");
                print_aligned(fp, "theta_ref");

                if (erotgFitPOT == rotg->eFittype)
                {
                    /* Output the set of angles around the reference angle */
                    for (i = 0; i < rotg->PotAngle_nstep; i++)
                    {
                        sprintf(buf, "E(%g)", erg->PotAngleFit->degangle[i]);
                        print_aligned(fp, buf);
                    }
                }
                else
                {
                    /* Output fit angle for each slab */
                    print_aligned_short(fp, "slab");
                    print_aligned_short(fp, "atoms");
                    print_aligned(fp, "theta_fit");
                    print_aligned_short(fp, "slab");
                    print_aligned_short(fp, "atoms");
                    print_aligned(fp, "theta_fit");
                    fprintf(fp, " ...");
                }
                fprintf(fp, "\n");
            }
        }
        fflush(fp);
    }

    return fp;
}


/* Open torque output file and write some information about it's structure.
 * Call on master only */
static FILE *open_torque_out(const char *fn, t_rot *rot, const output_env_t oenv)
{
    FILE      *fp;
<<<<<<< HEAD
    int       g;
=======
    int        g;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_rotgrp  *rotg;


    if (rot->enfrot->Flags & MD_APPENDFILES)
    {
<<<<<<< HEAD
        fp = gmx_fio_fopen(fn,"a");
    }
    else
    {
        fp = open_output_file(fn, rot->nstsout,"torques");

        for (g=0; g<rot->ngrp; g++)
=======
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = open_output_file(fn, rot->nstsout, "torques");

        for (g = 0; g < rot->ngrp; g++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            rotg = &rot->grp[g];
            if (ISFLEX(rotg))
            {
                fprintf(fp, "# Rotation group %d (%s), slab distance %f nm.\n", g, erotg_names[rotg->eType], rotg->slab_dist);
                fprintf(fp, "# The scalar tau is the torque (kJ/mol) in the direction of the rotation vector.\n");
                fprintf(fp, "# To obtain the vectorial torque, multiply tau with\n");
                fprintf(fp, "# rot_vec%d            %10.3e %10.3e %10.3e\n", g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
                fprintf(fp, "#\n");
            }
        }
        fprintf(fp, "# Legend for the following data columns: (tau=torque for that slab):\n");
        fprintf(fp, "#     ");
        print_aligned_short(fp, "t");
        print_aligned_short(fp, "grp");
        print_aligned_short(fp, "slab");
        print_aligned(fp, "tau");
        print_aligned_short(fp, "slab");
        print_aligned(fp, "tau");
        fprintf(fp, " ...\n");
        fflush(fp);
    }

    return fp;
}


static void swap_val(double* vec, int i, int j)
{
    double tmp = vec[j];
<<<<<<< HEAD
    
    
    vec[j]=vec[i];
    vec[i]=tmp;
=======


    vec[j] = vec[i];
    vec[i] = tmp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static void swap_col(double **mat, int i, int j)
{
    double tmp[3] = {mat[0][j], mat[1][j], mat[2][j]};
<<<<<<< HEAD
    
    
    mat[0][j]=mat[0][i];
    mat[1][j]=mat[1][i];
    mat[2][j]=mat[2][i];
    
    mat[0][i]=tmp[0];
    mat[1][i]=tmp[1];
    mat[2][i]=tmp[2];
} 
=======


    mat[0][j] = mat[0][i];
    mat[1][j] = mat[1][i];
    mat[2][j] = mat[2][i];

    mat[0][i] = tmp[0];
    mat[1][i] = tmp[1];
    mat[2][i] = tmp[2];
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/* Eigenvectors are stored in columns of eigen_vec */
static void diagonalize_symmetric(
        double **matrix,
        double **eigen_vec,
<<<<<<< HEAD
        double eigenval[3])
{
    int n_rot;
    
    
    jacobi(matrix,3,eigenval,eigen_vec,&n_rot);
    
=======
        double   eigenval[3])
{
    int n_rot;


    jacobi(matrix, 3, eigenval, eigen_vec, &n_rot);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* sort in ascending order */
    if (eigenval[0] > eigenval[1])
    {
        swap_val(eigenval, 0, 1);
        swap_col(eigen_vec, 0, 1);
<<<<<<< HEAD
    } 
=======
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (eigenval[1] > eigenval[2])
    {
        swap_val(eigenval, 1, 2);
        swap_col(eigen_vec, 1, 2);
    }
    if (eigenval[0] > eigenval[1])
    {
        swap_val(eigenval, 0, 1);
        swap_col(eigen_vec, 0, 1);
    }
}


static void align_with_z(
        rvec* s,           /* Structure to align */
<<<<<<< HEAD
        int natoms,
        rvec axis)
{
    int    i, j, k;
    rvec   zet = {0.0, 0.0, 1.0};
    rvec   rot_axis={0.0, 0.0, 0.0};
    rvec   *rotated_str=NULL;
    real   ooanorm;
    real   angle;
    matrix rotmat;
    
    
=======
        int   natoms,
        rvec  axis)
{
    int     i, j, k;
    rvec    zet         = {0.0, 0.0, 1.0};
    rvec    rot_axis    = {0.0, 0.0, 0.0};
    rvec   *rotated_str = NULL;
    real    ooanorm;
    real    angle;
    matrix  rotmat;


>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    snew(rotated_str, natoms);

    /* Normalize the axis */
    ooanorm = 1.0/norm(axis);
    svmul(ooanorm, axis, axis);
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Calculate the angle for the fitting procedure */
    cprod(axis, zet, rot_axis);
    angle = acos(axis[2]);
    if (angle < 0.0)
<<<<<<< HEAD
        angle += M_PI;
    
    /* Calculate the rotation matrix */
    calc_rotmat(rot_axis, angle*180.0/M_PI, rotmat);
    
    /* Apply the rotation matrix to s */
    for (i=0; i<natoms; i++)
    {    
        for(j=0; j<3; j++)
        {
            for(k=0; k<3; k++)
=======
    {
        angle += M_PI;
    }

    /* Calculate the rotation matrix */
    calc_rotmat(rot_axis, angle*180.0/M_PI, rotmat);

    /* Apply the rotation matrix to s */
    for (i = 0; i < natoms; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                rotated_str[i][j] += rotmat[j][k]*s[i][k];
            }
        }
    }
<<<<<<< HEAD
    
    /* Rewrite the rotated structure to s */
    for(i=0; i<natoms; i++)
    {
        for(j=0; j<3; j++)
        {
            s[i][j]=rotated_str[i][j];
        }
    }
    
    sfree(rotated_str);
} 


static void calc_correl_matrix(rvec* Xstr, rvec* Ystr, double** Rmat, int natoms)
{    
    int i, j, k;
 
    
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            Rmat[i][j] = 0.0;
    
    for (i=0; i<3; i++) 
        for (j=0; j<3; j++) 
            for (k=0; k<natoms; k++) 
                Rmat[i][j] += Ystr[k][i] * Xstr[k][j];
=======

    /* Rewrite the rotated structure to s */
    for (i = 0; i < natoms; i++)
    {
        for (j = 0; j < 3; j++)
        {
            s[i][j] = rotated_str[i][j];
        }
    }

    sfree(rotated_str);
}


static void calc_correl_matrix(rvec* Xstr, rvec* Ystr, double** Rmat, int natoms)
{
    int i, j, k;


    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            Rmat[i][j] = 0.0;
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < natoms; k++)
            {
                Rmat[i][j] += Ystr[k][i] * Xstr[k][j];
            }
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static void weigh_coords(rvec* str, real* weight, int natoms)
{
    int i, j;
<<<<<<< HEAD
    
    
    for(i=0; i<natoms; i++)
    {
        for(j=0; j<3; j++)
            str[i][j] *= sqrt(weight[i]);
    }  
=======


    for (i = 0; i < natoms; i++)
    {
        for (j = 0; j < 3; j++)
        {
            str[i][j] *= sqrt(weight[i]);
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static real opt_angle_analytic(
        rvec* ref_s,
        rvec* act_s,
<<<<<<< HEAD
        real* weight, 
        int natoms,
        rvec ref_com,
        rvec act_com,
        rvec axis)
{    
    int    i, j, k;
    rvec   *ref_s_1=NULL;
    rvec   *act_s_1=NULL;
    rvec   shift;
    double **Rmat, **RtR, **eigvec;
    double eigval[3];
    double V[3][3], WS[3][3];
    double rot_matrix[3][3];
    double opt_angle;
    
    
    /* Do not change the original coordinates */ 
    snew(ref_s_1, natoms);
    snew(act_s_1, natoms);
    for(i=0; i<natoms; i++)
=======
        real* weight,
        int   natoms,
        rvec  ref_com,
        rvec  act_com,
        rvec  axis)
{
    int      i, j, k;
    rvec    *ref_s_1 = NULL;
    rvec    *act_s_1 = NULL;
    rvec     shift;
    double **Rmat, **RtR, **eigvec;
    double   eigval[3];
    double   V[3][3], WS[3][3];
    double   rot_matrix[3][3];
    double   opt_angle;


    /* Do not change the original coordinates */
    snew(ref_s_1, natoms);
    snew(act_s_1, natoms);
    for (i = 0; i < natoms; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        copy_rvec(ref_s[i], ref_s_1[i]);
        copy_rvec(act_s[i], act_s_1[i]);
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Translate the structures to the origin */
    shift[XX] = -ref_com[XX];
    shift[YY] = -ref_com[YY];
    shift[ZZ] = -ref_com[ZZ];
    translate_x(ref_s_1, natoms, shift);
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    shift[XX] = -act_com[XX];
    shift[YY] = -act_com[YY];
    shift[ZZ] = -act_com[ZZ];
    translate_x(act_s_1, natoms, shift);
<<<<<<< HEAD
    
    /* Align rotation axis with z */
    align_with_z(ref_s_1, natoms, axis);
    align_with_z(act_s_1, natoms, axis);
    
    /* Correlation matrix */
    Rmat = allocate_square_matrix(3);
    
    for (i=0; i<natoms; i++)
    {
        ref_s_1[i][2]=0.0;
        act_s_1[i][2]=0.0;
    }
    
=======

    /* Align rotation axis with z */
    align_with_z(ref_s_1, natoms, axis);
    align_with_z(act_s_1, natoms, axis);

    /* Correlation matrix */
    Rmat = allocate_square_matrix(3);

    for (i = 0; i < natoms; i++)
    {
        ref_s_1[i][2] = 0.0;
        act_s_1[i][2] = 0.0;
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Weight positions with sqrt(weight) */
    if (NULL != weight)
    {
        weigh_coords(ref_s_1, weight, natoms);
        weigh_coords(act_s_1, weight, natoms);
    }
<<<<<<< HEAD
    
    /* Calculate correlation matrices R=YXt (X=ref_s; Y=act_s) */
    calc_correl_matrix(ref_s_1, act_s_1, Rmat, natoms);
    
    /* Calculate RtR */
    RtR = allocate_square_matrix(3);
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            for (k=0; k<3; k++)
=======

    /* Calculate correlation matrices R=YXt (X=ref_s; Y=act_s) */
    calc_correl_matrix(ref_s_1, act_s_1, Rmat, natoms);

    /* Calculate RtR */
    RtR = allocate_square_matrix(3);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                RtR[i][j] += Rmat[k][i] * Rmat[k][j];
            }
        }
    }
    /* Diagonalize RtR */
<<<<<<< HEAD
    snew(eigvec,3);
    for (i=0; i<3; i++)
        snew(eigvec[i],3);
    
    diagonalize_symmetric(RtR, eigvec, eigval);
    swap_col(eigvec,0,1);
    swap_col(eigvec,1,2);
    swap_val(eigval,0,1);
    swap_val(eigval,1,2);
    
    /* Calculate V */
    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
=======
    snew(eigvec, 3);
    for (i = 0; i < 3; i++)
    {
        snew(eigvec[i], 3);
    }

    diagonalize_symmetric(RtR, eigvec, eigval);
    swap_col(eigvec, 0, 1);
    swap_col(eigvec, 1, 2);
    swap_val(eigval, 0, 1);
    swap_val(eigval, 1, 2);

    /* Calculate V */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            V[i][j]  = 0.0;
            WS[i][j] = 0.0;
        }
    }
<<<<<<< HEAD
    
    for (i=0; i<2; i++)
        for (j=0; j<2; j++)
            WS[i][j] = eigvec[i][j] / sqrt(eigval[j]);
    
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            for (k=0; k<3; k++)
=======

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            WS[i][j] = eigvec[i][j] / sqrt(eigval[j]);
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                V[i][j] += Rmat[i][k]*WS[k][j];
            }
        }
    }
    free_square_matrix(Rmat, 3);
<<<<<<< HEAD
    
    /* Calculate optimal rotation matrix */
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            rot_matrix[i][j] = 0.0;
    
    for (i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            for(k=0; k<3; k++){
=======

    /* Calculate optimal rotation matrix */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            rot_matrix[i][j] = 0.0;
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                rot_matrix[i][j] += eigvec[i][k]*V[j][k];
            }
        }
    }
    rot_matrix[2][2] = 1.0;
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* In some cases abs(rot_matrix[0][0]) can be slighly larger
     * than unity due to numerical inacurracies. To be able to calculate
     * the acos function, we put these values back in range. */
    if (rot_matrix[0][0] > 1.0)
    {
        rot_matrix[0][0] = 1.0;
    }
    else if (rot_matrix[0][0] < -1.0)
    {
        rot_matrix[0][0] = -1.0;
    }

    /* Determine the optimal rotation angle: */
    opt_angle = (-1.0)*acos(rot_matrix[0][0])*180.0/M_PI;
    if (rot_matrix[0][1] < 0.0)
<<<<<<< HEAD
        opt_angle = (-1.0)*opt_angle;
        
=======
    {
        opt_angle = (-1.0)*opt_angle;
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Give back some memory */
    free_square_matrix(RtR, 3);
    sfree(ref_s_1);
    sfree(act_s_1);
<<<<<<< HEAD
    for (i=0; i<3; i++)
        sfree(eigvec[i]);
    sfree(eigvec);
    
=======
    for (i = 0; i < 3; i++)
    {
        sfree(eigvec[i]);
    }
    sfree(eigvec);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return (real) opt_angle;
}


/* Determine angle of the group by RMSD fit to the reference */
/* Not parallelized, call this routine only on the master */
static real flex_fit_angle(t_rotgrp *rotg)
{
<<<<<<< HEAD
    int         i;
    rvec        *fitcoords=NULL;
    rvec        center;         /* Center of positions passed to the fit routine */
    real        fitangle;       /* Angle of the rotation group derived by fitting */
    rvec        coord;
    real        scal;
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;
=======
    int             i;
    rvec           *fitcoords = NULL;
    rvec            center;     /* Center of positions passed to the fit routine */
    real            fitangle;   /* Angle of the rotation group derived by fitting */
    rvec            coord;
    real            scal;
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Get the center of the rotation group.
     * Note, again, erg->xc has been sorted in do_flexible */
    get_center(erg->xc, erg->mc_sorted, rotg->nat, center);

    /* === Determine the optimal fit angle for the rotation group === */
    if (rotg->eFittype == erotgFitNORM)
    {
        /* Normalize every position to it's reference length */
<<<<<<< HEAD
        for (i=0; i<rotg->nat; i++)
=======
        for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* Put the center of the positions into the origin */
            rvec_sub(erg->xc[i], center, coord);
            /* Determine the scaling factor for the length: */
            scal = erg->xc_ref_length[erg->xc_sortind[i]] / norm(coord);
            /* Get position, multiply with the scaling factor and save  */
            svmul(scal, coord, erg->xc_norm[i]);
        }
        fitcoords = erg->xc_norm;
    }
    else
    {
        fitcoords = erg->xc;
    }
    /* From the point of view of the current positions, the reference has rotated
     * backwards. Since we output the angle relative to the fixed reference,
     * we need the minus sign. */
    fitangle = -opt_angle_analytic(erg->xc_ref_sorted, fitcoords, erg->mc_sorted,
                                   rotg->nat, erg->xc_ref_center, center, rotg->vec);

    return fitangle;
}


/* Determine actual angle of each slab by RMSD fit to the reference */
/* Not parallelized, call this routine only on the master */
static void flex_fit_angle_perslab(
<<<<<<< HEAD
        int  g,
        t_rotgrp *rotg,
        double t,
        real degangle,
        FILE *fp)
{
    int         i,l,n,islab,ind;
    rvec        curr_x, ref_x;
    rvec        act_center;  /* Center of actual positions that are passed to the fit routine */
    rvec        ref_center;  /* Same for the reference positions */
    real        fitangle;    /* Angle of a slab derived from an RMSD fit to
                              * the reference structure at t=0  */
    t_gmx_slabdata *sd;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */
    real        OOm_av;      /* 1/average_mass of a rotation group atom */
    real        m_rel;       /* Relative mass of a rotation group atom  */


    erg=rotg->enfrotgrp;
=======
        int       g,
        t_rotgrp *rotg,
        double    t,
        real      degangle,
        FILE     *fp)
{
    int             i, l, n, islab, ind;
    rvec            curr_x, ref_x;
    rvec            act_center; /* Center of actual positions that are passed to the fit routine */
    rvec            ref_center; /* Same for the reference positions */
    real            fitangle;   /* Angle of a slab derived from an RMSD fit to
                                 * the reference structure at t=0  */
    t_gmx_slabdata *sd;
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data */
    real            OOm_av;     /* 1/average_mass of a rotation group atom */
    real            m_rel;      /* Relative mass of a rotation group atom  */


    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Average mass of a rotation group atom: */
    OOm_av = erg->invmass*rotg->nat;

    /**********************************/
    /* First collect the data we need */
    /**********************************/

    /* Collect the data for the individual slabs */
    for (n = erg->slab_first; n <= erg->slab_last; n++)
    {
<<<<<<< HEAD
        islab = n - erg->slab_first; /* slab index */
        sd = &(rotg->enfrotgrp->slab_data[islab]);
        sd->nat = erg->lastatom[islab]-erg->firstatom[islab]+1;
        ind = 0;

        /* Loop over the relevant atoms in the slab */
        for (l=erg->firstatom[islab]; l<=erg->lastatom[islab]; l++)
=======
        islab   = n - erg->slab_first; /* slab index */
        sd      = &(rotg->enfrotgrp->slab_data[islab]);
        sd->nat = erg->lastatom[islab]-erg->firstatom[islab]+1;
        ind     = 0;

        /* Loop over the relevant atoms in the slab */
        for (l = erg->firstatom[islab]; l <= erg->lastatom[islab]; l++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* Current position of this atom: x[ii][XX/YY/ZZ] */
            copy_rvec(erg->xc[l], curr_x);

            /* The (unrotated) reference position of this atom is copied to ref_x.
             * Beware, the xc coords have been sorted in do_flexible */
            copy_rvec(erg->xc_ref_sorted[l], ref_x);

            /* Save data for doing angular RMSD fit later */
            /* Save the current atom position */
            copy_rvec(curr_x, sd->x[ind]);
            /* Save the corresponding reference position */
<<<<<<< HEAD
            copy_rvec(ref_x , sd->ref[ind]);
=======
            copy_rvec(ref_x, sd->ref[ind]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            /* Maybe also mass-weighting was requested. If yes, additionally
             * multiply the weights with the relative mass of the atom. If not,
             * multiply with unity. */
            m_rel = erg->mc_sorted[l]*OOm_av;

            /* Save the weight for this atom in this slab */
            sd->weight[ind] = gaussian_weight(curr_x, rotg, n) * m_rel;

            /* Next atom in this slab */
            ind++;
        }
    }

    /******************************/
    /* Now do the fit calculation */
    /******************************/

    fprintf(fp, "%12.3e%6d%12.3f", t, g, degangle);

    /* === Now do RMSD fitting for each slab === */
    /* We require at least SLAB_MIN_ATOMS in a slab, such that the fit makes sense. */
#define SLAB_MIN_ATOMS 4

    for (n = erg->slab_first; n <= erg->slab_last; n++)
    {
        islab = n - erg->slab_first; /* slab index */
<<<<<<< HEAD
        sd = &(rotg->enfrotgrp->slab_data[islab]);
=======
        sd    = &(rotg->enfrotgrp->slab_data[islab]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (sd->nat >= SLAB_MIN_ATOMS)
        {
            /* Get the center of the slabs reference and current positions */
            get_center(sd->ref, sd->weight, sd->nat, ref_center);
<<<<<<< HEAD
            get_center(sd->x  , sd->weight, sd->nat, act_center);
=======
            get_center(sd->x, sd->weight, sd->nat, act_center);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (rotg->eFittype == erotgFitNORM)
            {
                /* Normalize every position to it's reference length
                 * prior to performing the fit */
<<<<<<< HEAD
                for (i=0; i<sd->nat;i++) /* Center */
                {
                    rvec_dec(sd->ref[i], ref_center);
                    rvec_dec(sd->x[i]  , act_center);
=======
                for (i = 0; i < sd->nat; i++) /* Center */
                {
                    rvec_dec(sd->ref[i], ref_center);
                    rvec_dec(sd->x[i], act_center);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    /* Normalize x_i such that it gets the same length as ref_i */
                    svmul( norm(sd->ref[i])/norm(sd->x[i]), sd->x[i], sd->x[i] );
                }
                /* We already subtracted the centers */
                clear_rvec(ref_center);
                clear_rvec(act_center);
            }
            fitangle = -opt_angle_analytic(sd->ref, sd->x, sd->weight, sd->nat,
                                           ref_center, act_center, rotg->vec);
            fprintf(fp, "%6d%6d%12.3f", n, sd->nat, fitangle);
        }
    }
<<<<<<< HEAD
    fprintf(fp     , "\n");
=======
    fprintf(fp, "\n");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#undef SLAB_MIN_ATOMS
}


/* Shift x with is */
static gmx_inline void shift_single_coord(matrix box, rvec x, const ivec is)
{
<<<<<<< HEAD
    int tx,ty,tz;


    tx=is[XX];
    ty=is[YY];
    tz=is[ZZ];

    if(TRICLINIC(box))
=======
    int tx, ty, tz;


    tx = is[XX];
    ty = is[YY];
    tz = is[ZZ];

    if (TRICLINIC(box))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        x[XX] += tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
        x[YY] += ty*box[YY][YY]+tz*box[ZZ][YY];
        x[ZZ] += tz*box[ZZ][ZZ];
<<<<<<< HEAD
    } else
=======
    }
    else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        x[XX] += tx*box[XX][XX];
        x[YY] += ty*box[YY][YY];
        x[ZZ] += tz*box[ZZ][ZZ];
    }
}


/* Determine the 'home' slab of this atom which is the
 * slab with the highest Gaussian weight of all */
#define round(a) (int)(a+0.5)
static gmx_inline int get_homeslab(
<<<<<<< HEAD
        rvec curr_x,   /* The position for which the home slab shall be determined */ 
=======
        rvec curr_x,   /* The position for which the home slab shall be determined */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        rvec rotvec,   /* The rotation vector */
        real slabdist) /* The slab distance */
{
    real dist;
<<<<<<< HEAD
    
    
    /* The distance of the atom to the coordinate center (where the
     * slab with index 0) is */
    dist = iprod(rotvec, curr_x);
    
    return round(dist / slabdist); 
=======


    /* The distance of the atom to the coordinate center (where the
     * slab with index 0) is */
    dist = iprod(rotvec, curr_x);

    return round(dist / slabdist);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* For a local atom determine the relevant slabs, i.e. slabs in
 * which the gaussian is larger than min_gaussian
 */
static int get_single_atom_gaussians(
<<<<<<< HEAD
        rvec      curr_x,
        t_rotgrp  *rotg)
{
   int slab, homeslab;
   real g;
   int count = 0;
   gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */

   
   erg=rotg->enfrotgrp;
   
   /* Determine the 'home' slab of this atom: */
   homeslab = get_homeslab(curr_x, rotg->vec, rotg->slab_dist);

   /* First determine the weight in the atoms home slab: */
   g = gaussian_weight(curr_x, rotg, homeslab);
   
   erg->gn_atom[count] = g;
   erg->gn_slabind[count] = homeslab;
   count++;
   
   
   /* Determine the max slab */
   slab = homeslab;
   while (g > rotg->min_gaussian)
   {
       slab++;
       g = gaussian_weight(curr_x, rotg, slab);
       erg->gn_slabind[count]=slab;
       erg->gn_atom[count]=g;
       count++;
   }
   count--;
   
   /* Determine the max slab */
   slab = homeslab;
   do
   {
       slab--;
       g = gaussian_weight(curr_x, rotg, slab);       
       erg->gn_slabind[count]=slab;
       erg->gn_atom[count]=g;
       count++;
   }
   while (g > rotg->min_gaussian);
   count--;
   
   return count;
=======
        rvec       curr_x,
        t_rotgrp  *rotg)
{
    int             slab, homeslab;
    real            g;
    int             count = 0;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;

    /* Determine the 'home' slab of this atom: */
    homeslab = get_homeslab(curr_x, rotg->vec, rotg->slab_dist);

    /* First determine the weight in the atoms home slab: */
    g = gaussian_weight(curr_x, rotg, homeslab);

    erg->gn_atom[count]    = g;
    erg->gn_slabind[count] = homeslab;
    count++;


    /* Determine the max slab */
    slab = homeslab;
    while (g > rotg->min_gaussian)
    {
        slab++;
        g = gaussian_weight(curr_x, rotg, slab);
        erg->gn_slabind[count] = slab;
        erg->gn_atom[count]    = g;
        count++;
    }
    count--;

    /* Determine the max slab */
    slab = homeslab;
    do
    {
        slab--;
        g = gaussian_weight(curr_x, rotg, slab);
        erg->gn_slabind[count] = slab;
        erg->gn_atom[count]    = g;
        count++;
    }
    while (g > rotg->min_gaussian);
    count--;

    return count;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static void flex2_precalc_inner_sum(t_rotgrp *rotg)
{
<<<<<<< HEAD
    int  i,n,islab;
    rvec  xi;                /* positions in the i-sum                        */
    rvec  xcn, ycn;          /* the current and the reference slab centers    */
    real gaussian_xi;
    rvec yi0;
    rvec  rin;               /* Helper variables                              */
    real  fac,fac2;
    rvec innersumvec;
    real OOpsii,OOpsiistar;
    real sin_rin;          /* s_ii.r_ii */
    rvec s_in,tmpvec,tmpvec2;
    real mi,wi;            /* Mass-weighting of the positions                 */
    real N_M;              /* N/M                                             */
    gmx_enfrotgrp_t erg;    /* Pointer to enforced rotation group data */


    erg=rotg->enfrotgrp;
    N_M = rotg->nat * erg->invmass;

    /* Loop over all slabs that contain something */
    for (n=erg->slab_first; n <= erg->slab_last; n++)
=======
    int             i, n, islab;
    rvec            xi;       /* positions in the i-sum                        */
    rvec            xcn, ycn; /* the current and the reference slab centers    */
    real            gaussian_xi;
    rvec            yi0;
    rvec            rin;     /* Helper variables                              */
    real            fac, fac2;
    rvec            innersumvec;
    real            OOpsii, OOpsiistar;
    real            sin_rin; /* s_ii.r_ii */
    rvec            s_in, tmpvec, tmpvec2;
    real            mi, wi;  /* Mass-weighting of the positions                 */
    real            N_M;     /* N/M                                             */
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;
    N_M = rotg->nat * erg->invmass;

    /* Loop over all slabs that contain something */
    for (n = erg->slab_first; n <= erg->slab_last; n++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        islab = n - erg->slab_first; /* slab index */

        /* The current center of this slab is saved in xcn: */
        copy_rvec(erg->slab_center[islab], xcn);
        /* ... and the reference center in ycn: */
        copy_rvec(erg->slab_center_ref[islab+erg->slab_buffer], ycn);

        /*** D. Calculate the whole inner sum used for second and third sum */
        /* For slab n, we need to loop over all atoms i again. Since we sorted
         * the atoms with respect to the rotation vector, we know that it is sufficient
         * to calculate from firstatom to lastatom only. All other contributions will
         * be very small. */
        clear_rvec(innersumvec);
        for (i = erg->firstatom[islab]; i <= erg->lastatom[islab]; i++)
        {
            /* Coordinate xi of this atom */
<<<<<<< HEAD
            copy_rvec(erg->xc[i],xi);

            /* The i-weights */
            gaussian_xi = gaussian_weight(xi,rotg,n);
            mi = erg->mc_sorted[i];  /* need the sorted mass here */
            wi = N_M*mi;

            /* Calculate rin */
            copy_rvec(erg->xc_ref_sorted[i],yi0); /* Reference position yi0   */
            rvec_sub(yi0, ycn, tmpvec2);          /* tmpvec2 = yi0 - ycn      */
            mvmul(erg->rotmat, tmpvec2, rin);     /* rin = Omega.(yi0 - ycn)  */
=======
            copy_rvec(erg->xc[i], xi);

            /* The i-weights */
            gaussian_xi = gaussian_weight(xi, rotg, n);
            mi          = erg->mc_sorted[i]; /* need the sorted mass here */
            wi          = N_M*mi;

            /* Calculate rin */
            copy_rvec(erg->xc_ref_sorted[i], yi0); /* Reference position yi0   */
            rvec_sub(yi0, ycn, tmpvec2);           /* tmpvec2 = yi0 - ycn      */
            mvmul(erg->rotmat, tmpvec2, rin);      /* rin = Omega.(yi0 - ycn)  */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            /* Calculate psi_i* and sin */
            rvec_sub(xi, xcn, tmpvec2);           /* tmpvec2 = xi - xcn       */
            cprod(rotg->vec, tmpvec2, tmpvec);    /* tmpvec = v x (xi - xcn)  */
            OOpsiistar = norm2(tmpvec)+rotg->eps; /* OOpsii* = 1/psii* = |v x (xi-xcn)|^2 + eps */
<<<<<<< HEAD
            OOpsii = norm(tmpvec);                /* OOpsii = 1 / psii = |v x (xi - xcn)| */

                                       /*         v x (xi - xcn)          */
            unitv(tmpvec, s_in);       /*  sin = ----------------         */
                                       /*        |v x (xi - xcn)|         */

            sin_rin=iprod(s_in,rin);   /* sin_rin = sin . rin             */
=======
            OOpsii     = norm(tmpvec);            /* OOpsii = 1 / psii = |v x (xi - xcn)| */

            /*                           *         v x (xi - xcn)          */
            unitv(tmpvec, s_in);        /*  sin = ----------------         */
                                        /*        |v x (xi - xcn)|         */

            sin_rin = iprod(s_in, rin); /* sin_rin = sin . rin             */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            /* Now the whole sum */
            fac = OOpsii/OOpsiistar;
            svmul(fac, rin, tmpvec);
            fac2 = fac*fac*OOpsii;
            svmul(fac2*sin_rin, s_in, tmpvec2);
            rvec_dec(tmpvec, tmpvec2);

            svmul(wi*gaussian_xi*sin_rin, tmpvec, tmpvec2);

<<<<<<< HEAD
            rvec_inc(innersumvec,tmpvec2);
=======
            rvec_inc(innersumvec, tmpvec2);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        } /* now we have the inner sum, used both for sum2 and sum3 */

        /* Save it to be used in do_flex2_lowlevel */
        copy_rvec(innersumvec, erg->slab_innersumvec[islab]);
    } /* END of loop over slabs */
}


static void flex_precalc_inner_sum(t_rotgrp *rotg)
{
<<<<<<< HEAD
    int   i,n,islab;
    rvec  xi;                /* position                                      */
    rvec  xcn, ycn;          /* the current and the reference slab centers    */
    rvec  qin,rin;           /* q_i^n and r_i^n                               */
    real  bin;
    rvec  tmpvec;
    rvec  innersumvec;       /* Inner part of sum_n2                          */
    real  gaussian_xi;       /* Gaussian weight gn(xi)                        */
    real  mi,wi;             /* Mass-weighting of the positions               */
    real  N_M;               /* N/M                                           */

    gmx_enfrotgrp_t erg;    /* Pointer to enforced rotation group data */


    erg=rotg->enfrotgrp;
    N_M = rotg->nat * erg->invmass;

    /* Loop over all slabs that contain something */
    for (n=erg->slab_first; n <= erg->slab_last; n++)
=======
    int             i, n, islab;
    rvec            xi;       /* position                                      */
    rvec            xcn, ycn; /* the current and the reference slab centers    */
    rvec            qin, rin; /* q_i^n and r_i^n                               */
    real            bin;
    rvec            tmpvec;
    rvec            innersumvec; /* Inner part of sum_n2                          */
    real            gaussian_xi; /* Gaussian weight gn(xi)                        */
    real            mi, wi;      /* Mass-weighting of the positions               */
    real            N_M;         /* N/M                                           */

    gmx_enfrotgrp_t erg;         /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;
    N_M = rotg->nat * erg->invmass;

    /* Loop over all slabs that contain something */
    for (n = erg->slab_first; n <= erg->slab_last; n++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        islab = n - erg->slab_first; /* slab index */

        /* The current center of this slab is saved in xcn: */
        copy_rvec(erg->slab_center[islab], xcn);
        /* ... and the reference center in ycn: */
        copy_rvec(erg->slab_center_ref[islab+erg->slab_buffer], ycn);

        /* For slab n, we need to loop over all atoms i again. Since we sorted
         * the atoms with respect to the rotation vector, we know that it is sufficient
         * to calculate from firstatom to lastatom only. All other contributions will
         * be very small. */
        clear_rvec(innersumvec);
<<<<<<< HEAD
        for (i=erg->firstatom[islab]; i<=erg->lastatom[islab]; i++)
        {
            /* Coordinate xi of this atom */
            copy_rvec(erg->xc[i],xi);

            /* The i-weights */
            gaussian_xi = gaussian_weight(xi,rotg,n);
            mi = erg->mc_sorted[i];  /* need the sorted mass here */
            wi = N_M*mi;

            /* Calculate rin and qin */
            rvec_sub(erg->xc_ref_sorted[i], ycn, tmpvec); /* tmpvec = yi0-ycn */
            mvmul(erg->rotmat, tmpvec, rin);      /* rin = Omega.(yi0 - ycn)  */
            cprod(rotg->vec, rin, tmpvec);    /* tmpvec = v x Omega*(yi0-ycn) */

                                             /*        v x Omega*(yi0-ycn)    */
=======
        for (i = erg->firstatom[islab]; i <= erg->lastatom[islab]; i++)
        {
            /* Coordinate xi of this atom */
            copy_rvec(erg->xc[i], xi);

            /* The i-weights */
            gaussian_xi = gaussian_weight(xi, rotg, n);
            mi          = erg->mc_sorted[i]; /* need the sorted mass here */
            wi          = N_M*mi;

            /* Calculate rin and qin */
            rvec_sub(erg->xc_ref_sorted[i], ycn, tmpvec); /* tmpvec = yi0-ycn */
            mvmul(erg->rotmat, tmpvec, rin);              /* rin = Omega.(yi0 - ycn)  */
            cprod(rotg->vec, rin, tmpvec);                /* tmpvec = v x Omega*(yi0-ycn) */

            /*                                *        v x Omega*(yi0-ycn)    */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            unitv(tmpvec, qin);              /* qin = ---------------------   */
                                             /*       |v x Omega*(yi0-ycn)|   */

            /* Calculate bin */
            rvec_sub(xi, xcn, tmpvec);            /* tmpvec = xi-xcn          */
            bin = iprod(qin, tmpvec);             /* bin  = qin*(xi-xcn)      */

            svmul(wi*gaussian_xi*bin, qin, tmpvec);

            /* Add this contribution to the inner sum: */
            rvec_add(innersumvec, tmpvec, innersumvec);
        } /* now we have the inner sum vector S^n for this slab */
<<<<<<< HEAD
        /* Save it to be used in do_flex_lowlevel */
=======
          /* Save it to be used in do_flex_lowlevel */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        copy_rvec(innersumvec, erg->slab_innersumvec[islab]);
    }
}


static real do_flex2_lowlevel(
        t_rotgrp  *rotg,
<<<<<<< HEAD
        real      sigma,    /* The Gaussian width sigma */
        rvec      x[],
        gmx_bool  bOutstepRot,
        gmx_bool  bOutstepSlab,
        matrix    box)
{
    int  count,ic,ii,j,m,n,islab,iigrp,ifit;
    rvec xj;                 /* position in the i-sum                         */
    rvec yj0;                /* the reference position in the j-sum           */
    rvec xcn, ycn;           /* the current and the reference slab centers    */
    real V;                  /* This node's part of the rotation pot. energy  */
    real gaussian_xj;        /* Gaussian weight                               */
    real beta;

    real  numerator,fit_numerator;
    rvec  rjn,fit_rjn;       /* Helper variables                              */
    real  fac,fac2;

    real OOpsij,OOpsijstar;
    real OOsigma2;           /* 1/(sigma^2)                                   */
    real sjn_rjn;
    real betasigpsi;
    rvec sjn,tmpvec,tmpvec2,yj0_ycn;
    rvec sum1vec_part,sum1vec,sum2vec_part,sum2vec,sum3vec,sum4vec,innersumvec;
    real sum3,sum4;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data       */
    real mj,wj;              /* Mass-weighting of the positions               */
    real N_M;                /* N/M                                           */
    real Wjn;                /* g_n(x_j) m_j / Mjn                            */
    gmx_bool bCalcPotFit;
=======
        real       sigma,   /* The Gaussian width sigma */
        rvec       x[],
        gmx_bool   bOutstepRot,
        gmx_bool   bOutstepSlab,
        matrix     box)
{
    int             count, ic, ii, j, m, n, islab, iigrp, ifit;
    rvec            xj;          /* position in the i-sum                         */
    rvec            yj0;         /* the reference position in the j-sum           */
    rvec            xcn, ycn;    /* the current and the reference slab centers    */
    real            V;           /* This node's part of the rotation pot. energy  */
    real            gaussian_xj; /* Gaussian weight                               */
    real            beta;

    real            numerator, fit_numerator;
    rvec            rjn, fit_rjn; /* Helper variables                              */
    real            fac, fac2;

    real            OOpsij, OOpsijstar;
    real            OOsigma2; /* 1/(sigma^2)                                   */
    real            sjn_rjn;
    real            betasigpsi;
    rvec            sjn, tmpvec, tmpvec2, yj0_ycn;
    rvec            sum1vec_part, sum1vec, sum2vec_part, sum2vec, sum3vec, sum4vec, innersumvec;
    real            sum3, sum4;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data       */
    real            mj, wj;  /* Mass-weighting of the positions               */
    real            N_M;     /* N/M                                           */
    real            Wjn;     /* g_n(x_j) m_j / Mjn                            */
    gmx_bool        bCalcPotFit;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* To calculate the torque per slab */
    rvec slab_force;         /* Single force from slab n on one atom          */
    rvec slab_sum1vec_part;
<<<<<<< HEAD
    real slab_sum3part,slab_sum4part;
    rvec slab_sum1vec, slab_sum2vec, slab_sum3vec, slab_sum4vec;


    erg=rotg->enfrotgrp;
=======
    real slab_sum3part, slab_sum4part;
    rvec slab_sum1vec, slab_sum2vec, slab_sum3vec, slab_sum4vec;


    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Pre-calculate the inner sums, so that we do not have to calculate
     * them again for every atom */
    flex2_precalc_inner_sum(rotg);

<<<<<<< HEAD
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT==rotg->eFittype);
=======
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT == rotg->eFittype);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
<<<<<<< HEAD
    N_M = rotg->nat * erg->invmass;
    V = 0.0;
    OOsigma2 = 1.0 / (sigma*sigma);
    for (j=0; j<erg->nat_loc; j++)
=======
    N_M      = rotg->nat * erg->invmass;
    V        = 0.0;
    OOsigma2 = 1.0 / (sigma*sigma);
    for (j = 0; j < erg->nat_loc; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Local index of a rotation group atom  */
        ii = erg->ind_loc[j];
        /* Position of this atom in the collective array */
        iigrp = erg->xc_ref_ind[j];
        /* Mass-weighting */
        mj = erg->mc[iigrp];  /* need the unsorted mass here */
        wj = N_M*mj;
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Current position of this atom: x[ii][XX/YY/ZZ]
         * Note that erg->xc_center contains the center of mass in case the flex2-t
         * potential was chosen. For the flex2 potential erg->xc_center must be
         * zero. */
        rvec_sub(x[ii], erg->xc_center, xj);

        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

        /* Determine the slabs to loop over, i.e. the ones with contributions
         * larger than min_gaussian */
        count = get_single_atom_gaussians(xj, rotg);
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        clear_rvec(sum1vec_part);
        clear_rvec(sum2vec_part);
        sum3 = 0.0;
        sum4 = 0.0;
        /* Loop over the relevant slabs for this atom */
<<<<<<< HEAD
        for (ic=0; ic < count; ic++)  
        {
            n = erg->gn_slabind[ic];
            
=======
        for (ic = 0; ic < count; ic++)
        {
            n = erg->gn_slabind[ic];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Get the precomputed Gaussian value of curr_slab for curr_x */
            gaussian_xj = erg->gn_atom[ic];

            islab = n - erg->slab_first; /* slab index */
<<<<<<< HEAD
            
            /* The (unrotated) reference position of this atom is copied to yj0: */
            copy_rvec(rotg->x_ref[iigrp], yj0);

            beta = calc_beta(xj, rotg,n);
=======

            /* The (unrotated) reference position of this atom is copied to yj0: */
            copy_rvec(rotg->x_ref[iigrp], yj0);

            beta = calc_beta(xj, rotg, n);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            /* The current center of this slab is saved in xcn: */
            copy_rvec(erg->slab_center[islab], xcn);
            /* ... and the reference center in ycn: */
            copy_rvec(erg->slab_center_ref[islab+erg->slab_buffer], ycn);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            rvec_sub(yj0, ycn, yj0_ycn);          /* yj0_ycn = yj0 - ycn      */

            /* Rotate: */
            mvmul(erg->rotmat, yj0_ycn, rjn);     /* rjn = Omega.(yj0 - ycn)  */
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Subtract the slab center from xj */
            rvec_sub(xj, xcn, tmpvec2);           /* tmpvec2 = xj - xcn       */

            /* Calculate sjn */
            cprod(rotg->vec, tmpvec2, tmpvec);    /* tmpvec = v x (xj - xcn)  */

            OOpsijstar = norm2(tmpvec)+rotg->eps; /* OOpsij* = 1/psij* = |v x (xj-xcn)|^2 + eps */

            numerator = sqr(iprod(tmpvec, rjn));
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /*********************************/
            /* Add to the rotation potential */
            /*********************************/
            V += 0.5*rotg->k*wj*gaussian_xj*numerator/OOpsijstar;

            /* If requested, also calculate the potential for a set of angles
             * near the current reference angle */
            if (bCalcPotFit)
            {
                for (ifit = 0; ifit < rotg->PotAngle_nstep; ifit++)
                {
                    mvmul(erg->PotAngleFit->rotmat[ifit], yj0_ycn, fit_rjn);
<<<<<<< HEAD
                    fit_numerator = sqr(iprod(tmpvec, fit_rjn));
=======
                    fit_numerator              = sqr(iprod(tmpvec, fit_rjn));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    erg->PotAngleFit->V[ifit] += 0.5*rotg->k*wj*gaussian_xj*fit_numerator/OOpsijstar;
                }
            }

            /*************************************/
            /* Now calculate the force on atom j */
            /*************************************/

            OOpsij = norm(tmpvec);    /* OOpsij = 1 / psij = |v x (xj - xcn)| */

<<<<<<< HEAD
                                           /*         v x (xj - xcn)          */
            unitv(tmpvec, sjn);            /*  sjn = ----------------         */
                                           /*        |v x (xj - xcn)|         */

            sjn_rjn=iprod(sjn,rjn);        /* sjn_rjn = sjn . rjn             */
=======
            /*                              *         v x (xj - xcn)          */
            unitv(tmpvec, sjn);            /*  sjn = ----------------         */
                                           /*        |v x (xj - xcn)|         */

            sjn_rjn = iprod(sjn, rjn);     /* sjn_rjn = sjn . rjn             */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


            /*** A. Calculate the first of the four sum terms: ****************/
            fac = OOpsij/OOpsijstar;
            svmul(fac, rjn, tmpvec);
            fac2 = fac*fac*OOpsij;
            svmul(fac2*sjn_rjn, sjn, tmpvec2);
            rvec_dec(tmpvec, tmpvec2);
            fac2 = wj*gaussian_xj; /* also needed for sum4 */
            svmul(fac2*sjn_rjn, tmpvec, slab_sum1vec_part);
            /********************/
            /*** Add to sum1: ***/
            /********************/
            rvec_inc(sum1vec_part, slab_sum1vec_part); /* sum1 still needs to vector multiplied with v */

            /*** B. Calculate the forth of the four sum terms: ****************/
            betasigpsi = beta*OOsigma2*OOpsij; /* this is also needed for sum3 */
            /********************/
            /*** Add to sum4: ***/
            /********************/
            slab_sum4part = fac2*betasigpsi*fac*sjn_rjn*sjn_rjn; /* Note that fac is still valid from above */
<<<<<<< HEAD
            sum4 += slab_sum4part;
=======
            sum4         += slab_sum4part;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            /*** C. Calculate Wjn for second and third sum */
            /* Note that we can safely divide by slab_weights since we check in
             * get_slab_centers that it is non-zero. */
            Wjn = gaussian_xj*mj/erg->slab_weights[islab];

            /* We already have precalculated the inner sum for slab n */
            copy_rvec(erg->slab_innersumvec[islab], innersumvec);

            /* Weigh the inner sum vector with Wjn */
            svmul(Wjn, innersumvec, innersumvec);

            /*** E. Calculate the second of the four sum terms: */
            /********************/
            /*** Add to sum2: ***/
            /********************/
            rvec_inc(sum2vec_part, innersumvec); /* sum2 still needs to be vector crossproduct'ed with v */
<<<<<<< HEAD
            
            /*** F. Calculate the third of the four sum terms: */
            slab_sum3part = betasigpsi * iprod(sjn, innersumvec);
            sum3 += slab_sum3part; /* still needs to be multiplied with v */
=======

            /*** F. Calculate the third of the four sum terms: */
            slab_sum3part = betasigpsi * iprod(sjn, innersumvec);
            sum3         += slab_sum3part; /* still needs to be multiplied with v */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            /*** G. Calculate the torque on the local slab's axis: */
            if (bOutstepRot)
            {
                /* Sum1 */
                cprod(slab_sum1vec_part, rotg->vec, slab_sum1vec);
                /* Sum2 */
                cprod(innersumvec, rotg->vec, slab_sum2vec);
                /* Sum3 */
                svmul(slab_sum3part, rotg->vec, slab_sum3vec);
                /* Sum4 */
                svmul(slab_sum4part, rotg->vec, slab_sum4vec);

                /* The force on atom ii from slab n only: */
<<<<<<< HEAD
                for (m=0; m<DIM; m++)
                    slab_force[m] = rotg->k * (-slab_sum1vec[m] + slab_sum2vec[m] - slab_sum3vec[m] + 0.5*slab_sum4vec[m]);
=======
                for (m = 0; m < DIM; m++)
                {
                    slab_force[m] = rotg->k * (-slab_sum1vec[m] + slab_sum2vec[m] - slab_sum3vec[m] + 0.5*slab_sum4vec[m]);
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

                erg->slab_torque_v[islab] += torque(rotg->vec, slab_force, xj, xcn);
            }
        } /* END of loop over slabs */

        /* Construct the four individual parts of the vector sum: */
        cprod(sum1vec_part, rotg->vec, sum1vec);      /* sum1vec =   { } x v  */
        cprod(sum2vec_part, rotg->vec, sum2vec);      /* sum2vec =   { } x v  */
        svmul(sum3, rotg->vec, sum3vec);              /* sum3vec =   { } . v  */
        svmul(sum4, rotg->vec, sum4vec);              /* sum4vec =   { } . v  */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
<<<<<<< HEAD
        for (m=0; m<DIM; m++)
            erg->f_rot_loc[j][m] = rotg->k * (-sum1vec[m] + sum2vec[m] - sum3vec[m] + 0.5*sum4vec[m]);
=======
        for (m = 0; m < DIM; m++)
        {
            erg->f_rot_loc[j][m] = rotg->k * (-sum1vec[m] + sum2vec[m] - sum3vec[m] + 0.5*sum4vec[m]);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef SUM_PARTS
        fprintf(stderr, "sum1: %15.8f %15.8f %15.8f\n",    -rotg->k*sum1vec[XX],    -rotg->k*sum1vec[YY],    -rotg->k*sum1vec[ZZ]);
        fprintf(stderr, "sum2: %15.8f %15.8f %15.8f\n",     rotg->k*sum2vec[XX],     rotg->k*sum2vec[YY],     rotg->k*sum2vec[ZZ]);
        fprintf(stderr, "sum3: %15.8f %15.8f %15.8f\n",    -rotg->k*sum3vec[XX],    -rotg->k*sum3vec[YY],    -rotg->k*sum3vec[ZZ]);
        fprintf(stderr, "sum4: %15.8f %15.8f %15.8f\n", 0.5*rotg->k*sum4vec[XX], 0.5*rotg->k*sum4vec[YY], 0.5*rotg->k*sum4vec[ZZ]);
#endif

        PRINT_FORCE_J

    } /* END of loop over local atoms */

    return V;
}


static real do_flex_lowlevel(
        t_rotgrp *rotg,
        real      sigma,     /* The Gaussian width sigma                      */
        rvec      x[],
        gmx_bool  bOutstepRot,
        gmx_bool  bOutstepSlab,
        matrix    box)
{
<<<<<<< HEAD
    int   count,ic,ifit,ii,j,m,n,islab,iigrp;
    rvec  xj,yj0;            /* current and reference position                */
    rvec  xcn, ycn;          /* the current and the reference slab centers    */
    rvec  yj0_ycn;           /* yj0 - ycn                                     */
    rvec  xj_xcn;            /* xj - xcn                                      */
    rvec  qjn,fit_qjn;       /* q_i^n                                         */
    rvec  sum_n1,sum_n2;     /* Two contributions to the rotation force       */
    rvec  innersumvec;       /* Inner part of sum_n2                          */
    rvec  s_n;
    rvec  force_n;           /* Single force from slab n on one atom          */
    rvec  force_n1,force_n2; /* First and second part of force_n              */
    rvec  tmpvec,tmpvec2,tmp_f;   /* Helper variables                         */
    real  V;                 /* The rotation potential energy                 */
    real  OOsigma2;          /* 1/(sigma^2)                                   */
    real  beta;              /* beta_n(xj)                                    */
    real  bjn, fit_bjn;      /* b_j^n                                         */
    real  gaussian_xj;       /* Gaussian weight gn(xj)                        */
    real  betan_xj_sigma2;
    real  mj,wj;             /* Mass-weighting of the positions               */
    real  N_M;               /* N/M                                           */
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data       */
    gmx_bool bCalcPotFit;

    
    erg=rotg->enfrotgrp;
=======
    int             count, ic, ifit, ii, j, m, n, islab, iigrp;
    rvec            xj, yj0;                /* current and reference position                */
    rvec            xcn, ycn;               /* the current and the reference slab centers    */
    rvec            yj0_ycn;                /* yj0 - ycn                                     */
    rvec            xj_xcn;                 /* xj - xcn                                      */
    rvec            qjn, fit_qjn;           /* q_i^n                                         */
    rvec            sum_n1, sum_n2;         /* Two contributions to the rotation force       */
    rvec            innersumvec;            /* Inner part of sum_n2                          */
    rvec            s_n;
    rvec            force_n;                /* Single force from slab n on one atom          */
    rvec            force_n1, force_n2;     /* First and second part of force_n              */
    rvec            tmpvec, tmpvec2, tmp_f; /* Helper variables                              */
    real            V;                      /* The rotation potential energy                 */
    real            OOsigma2;               /* 1/(sigma^2)                                   */
    real            beta;                   /* beta_n(xj)                                    */
    real            bjn, fit_bjn;           /* b_j^n                                         */
    real            gaussian_xj;            /* Gaussian weight gn(xj)                        */
    real            betan_xj_sigma2;
    real            mj, wj;                 /* Mass-weighting of the positions               */
    real            N_M;                    /* N/M                                           */
    gmx_enfrotgrp_t erg;                    /* Pointer to enforced rotation group data       */
    gmx_bool        bCalcPotFit;


    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Pre-calculate the inner sums, so that we do not have to calculate
     * them again for every atom */
    flex_precalc_inner_sum(rotg);

<<<<<<< HEAD
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT==rotg->eFittype);
=======
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT == rotg->eFittype);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    OOsigma2 = 1.0/(sigma*sigma);
<<<<<<< HEAD
    N_M = rotg->nat * erg->invmass;
    V = 0.0;
    for (j=0; j<erg->nat_loc; j++)
=======
    N_M      = rotg->nat * erg->invmass;
    V        = 0.0;
    for (j = 0; j < erg->nat_loc; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Local index of a rotation group atom  */
        ii = erg->ind_loc[j];
        /* Position of this atom in the collective array */
        iigrp = erg->xc_ref_ind[j];
        /* Mass-weighting */
        mj = erg->mc[iigrp];  /* need the unsorted mass here */
        wj = N_M*mj;
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Current position of this atom: x[ii][XX/YY/ZZ]
         * Note that erg->xc_center contains the center of mass in case the flex-t
         * potential was chosen. For the flex potential erg->xc_center must be
         * zero. */
        rvec_sub(x[ii], erg->xc_center, xj);
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

        /* Determine the slabs to loop over, i.e. the ones with contributions
         * larger than min_gaussian */
        count = get_single_atom_gaussians(xj, rotg);

        clear_rvec(sum_n1);
        clear_rvec(sum_n2);

        /* Loop over the relevant slabs for this atom */
<<<<<<< HEAD
        for (ic=0; ic < count; ic++)  
        {
            n = erg->gn_slabind[ic];
                
=======
        for (ic = 0; ic < count; ic++)
        {
            n = erg->gn_slabind[ic];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Get the precomputed Gaussian for xj in slab n */
            gaussian_xj = erg->gn_atom[ic];

            islab = n - erg->slab_first; /* slab index */
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* The (unrotated) reference position of this atom is saved in yj0: */
            copy_rvec(rotg->x_ref[iigrp], yj0);

            beta = calc_beta(xj, rotg, n);

            /* The current center of this slab is saved in xcn: */
            copy_rvec(erg->slab_center[islab], xcn);
            /* ... and the reference center in ycn: */
            copy_rvec(erg->slab_center_ref[islab+erg->slab_buffer], ycn);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            rvec_sub(yj0, ycn, yj0_ycn); /* yj0_ycn = yj0 - ycn */

            /* Rotate: */
            mvmul(erg->rotmat, yj0_ycn, tmpvec2); /* tmpvec2= Omega.(yj0-ycn) */
<<<<<<< HEAD
            
            /* Subtract the slab center from xj */
            rvec_sub(xj, xcn, xj_xcn);           /* xj_xcn = xj - xcn         */
            
            /* Calculate qjn */
            cprod(rotg->vec, tmpvec2, tmpvec); /* tmpvec= v x Omega.(yj0-ycn) */

                                 /*         v x Omega.(yj0-ycn)    */
            unitv(tmpvec,qjn);   /*  qjn = ---------------------   */
                                 /*        |v x Omega.(yj0-ycn)|   */

            bjn = iprod(qjn, xj_xcn);   /* bjn = qjn * (xj - xcn) */
            
=======

            /* Subtract the slab center from xj */
            rvec_sub(xj, xcn, xj_xcn);           /* xj_xcn = xj - xcn         */

            /* Calculate qjn */
            cprod(rotg->vec, tmpvec2, tmpvec); /* tmpvec= v x Omega.(yj0-ycn) */

            /*                         *         v x Omega.(yj0-ycn)    */
            unitv(tmpvec, qjn);       /*  qjn = ---------------------   */
                                      /*        |v x Omega.(yj0-ycn)|   */

            bjn = iprod(qjn, xj_xcn); /* bjn = qjn * (xj - xcn) */

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /*********************************/
            /* Add to the rotation potential */
            /*********************************/
            V += 0.5*rotg->k*wj*gaussian_xj*sqr(bjn);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* If requested, also calculate the potential for a set of angles
             * near the current reference angle */
            if (bCalcPotFit)
            {
                for (ifit = 0; ifit < rotg->PotAngle_nstep; ifit++)
                {
                    /* As above calculate Omega.(yj0-ycn), now for the other angles */
                    mvmul(erg->PotAngleFit->rotmat[ifit], yj0_ycn, tmpvec2); /* tmpvec2= Omega.(yj0-ycn) */
                    /* As above calculate qjn */
<<<<<<< HEAD
                    cprod(rotg->vec, tmpvec2, tmpvec); /* tmpvec= v x Omega.(yj0-ycn) */
                                             /*             v x Omega.(yj0-ycn)    */
                    unitv(tmpvec,fit_qjn);   /*  fit_qjn = ---------------------   */
                                             /*            |v x Omega.(yj0-ycn)|   */
                    fit_bjn = iprod(fit_qjn, xj_xcn);   /* fit_bjn = fit_qjn * (xj - xcn) */
=======
                    cprod(rotg->vec, tmpvec2, tmpvec);                       /* tmpvec= v x Omega.(yj0-ycn) */
                    /*                                                        *             v x Omega.(yj0-ycn)    */
                    unitv(tmpvec, fit_qjn);                                  /*  fit_qjn = ---------------------   */
                                                                             /*            |v x Omega.(yj0-ycn)|   */
                    fit_bjn = iprod(fit_qjn, xj_xcn);                        /* fit_bjn = fit_qjn * (xj - xcn) */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    /* Add to the rotation potential for this angle */
                    erg->PotAngleFit->V[ifit] += 0.5*rotg->k*wj*gaussian_xj*sqr(fit_bjn);
                }
            }

            /****************************************************************/
            /* sum_n1 will typically be the main contribution to the force: */
            /****************************************************************/
            betan_xj_sigma2 = beta*OOsigma2;  /*  beta_n(xj)/sigma^2  */

            /* The next lines calculate
             *  qjn - (bjn*beta(xj)/(2sigma^2))v  */
            svmul(bjn*0.5*betan_xj_sigma2, rotg->vec, tmpvec2);
<<<<<<< HEAD
            rvec_sub(qjn,tmpvec2,tmpvec);

            /* Multiply with gn(xj)*bjn: */
            svmul(gaussian_xj*bjn,tmpvec,tmpvec2);

            /* Sum over n: */
            rvec_inc(sum_n1,tmpvec2);
            
            /* We already have precalculated the Sn term for slab n */
            copy_rvec(erg->slab_innersumvec[islab], s_n);
                                                                          /*          beta_n(xj)              */
=======
            rvec_sub(qjn, tmpvec2, tmpvec);

            /* Multiply with gn(xj)*bjn: */
            svmul(gaussian_xj*bjn, tmpvec, tmpvec2);

            /* Sum over n: */
            rvec_inc(sum_n1, tmpvec2);

            /* We already have precalculated the Sn term for slab n */
            copy_rvec(erg->slab_innersumvec[islab], s_n);
            /*                                                             *          beta_n(xj)              */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            svmul(betan_xj_sigma2*iprod(s_n, xj_xcn), rotg->vec, tmpvec); /* tmpvec = ---------- s_n (xj-xcn) */
                                                                          /*            sigma^2               */

            rvec_sub(s_n, tmpvec, innersumvec);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* We can safely divide by slab_weights since we check in get_slab_centers
             * that it is non-zero. */
            svmul(gaussian_xj/erg->slab_weights[islab], innersumvec, innersumvec);

            rvec_add(sum_n2, innersumvec, sum_n2);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Calculate the torque: */
            if (bOutstepRot)
            {
                /* The force on atom ii from slab n only: */
<<<<<<< HEAD
                svmul(-rotg->k*wj, tmpvec2    , force_n1); /* part 1 */
=======
                svmul(-rotg->k*wj, tmpvec2, force_n1);     /* part 1 */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                svmul( rotg->k*mj, innersumvec, force_n2); /* part 2 */
                rvec_add(force_n1, force_n2, force_n);
                erg->slab_torque_v[islab] += torque(rotg->vec, force_n, xj, xcn);
            }
        } /* END of loop over slabs */

        /* Put both contributions together: */
        svmul(wj, sum_n1, sum_n1);
        svmul(mj, sum_n2, sum_n2);
<<<<<<< HEAD
        rvec_sub(sum_n2,sum_n1,tmp_f); /* F = -grad V */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for(m=0; m<DIM; m++)
            erg->f_rot_loc[j][m] = rotg->k*tmp_f[m];
=======
        rvec_sub(sum_n2, sum_n1, tmp_f); /* F = -grad V */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (m = 0; m < DIM; m++)
        {
            erg->f_rot_loc[j][m] = rotg->k*tmp_f[m];
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        PRINT_FORCE_J

    } /* END of loop over local atoms */

    return V;
}

#ifdef PRINT_COORDS
static void print_coordinates(t_rotgrp *rotg, rvec x[], matrix box, int step)
{
<<<<<<< HEAD
    int i;
    static FILE *fp;
    static char buf[STRLEN];
    static gmx_bool bFirst=1;
=======
    int             i;
    static FILE    *fp;
    static char     buf[STRLEN];
    static gmx_bool bFirst = 1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    if (bFirst)
    {
        sprintf(buf, "coords%d.txt", cr->nodeid);
<<<<<<< HEAD
        fp = fopen(buf, "w");
=======
        fp     = fopen(buf, "w");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        bFirst = 0;
    }

    fprintf(fp, "\nStep %d\n", step);
    fprintf(fp, "box: %f %f %f %f %f %f %f %f %f\n",
            box[XX][XX], box[XX][YY], box[XX][ZZ],
            box[YY][XX], box[YY][YY], box[YY][ZZ],
            box[ZZ][XX], box[ZZ][ZZ], box[ZZ][ZZ]);
<<<<<<< HEAD
    for (i=0; i<rotg->nat; i++)
=======
    for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        fprintf(fp, "%4d  %f %f %f\n", i,
                erg->xc[i][XX], erg->xc[i][YY], erg->xc[i][ZZ]);
    }
    fflush(fp);

}
#endif


static int projection_compare(const void *a, const void *b)
{
    sort_along_vec_t *xca, *xcb;
<<<<<<< HEAD
    
    
    xca = (sort_along_vec_t *)a;
    xcb = (sort_along_vec_t *)b;
    
    if (xca->xcproj < xcb->xcproj)
        return -1;
    else if (xca->xcproj > xcb->xcproj)
        return 1;
    else
        return 0;
=======


    xca = (sort_along_vec_t *)a;
    xcb = (sort_along_vec_t *)b;

    if (xca->xcproj < xcb->xcproj)
    {
        return -1;
    }
    else if (xca->xcproj > xcb->xcproj)
    {
        return 1;
    }
    else
    {
        return 0;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static void sort_collective_coordinates(
<<<<<<< HEAD
        t_rotgrp *rotg,         /* Rotation group */
        sort_along_vec_t *data) /* Buffer for sorting the positions */
{
    int i;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;
    
    /* The projection of the position vector on the rotation vector is
     * the relevant value for sorting. Fill the 'data' structure */
    for (i=0; i<rotg->nat; i++)
=======
        t_rotgrp         *rotg, /* Rotation group */
        sort_along_vec_t *data) /* Buffer for sorting the positions */
{
    int             i;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;

    /* The projection of the position vector on the rotation vector is
     * the relevant value for sorting. Fill the 'data' structure */
    for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        data[i].xcproj = iprod(erg->xc[i], rotg->vec);  /* sort criterium */
        data[i].m      = erg->mc[i];
        data[i].ind    = i;
<<<<<<< HEAD
        copy_rvec(erg->xc[i]    , data[i].x    );
=======
        copy_rvec(erg->xc[i], data[i].x    );
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        copy_rvec(rotg->x_ref[i], data[i].x_ref);
    }
    /* Sort the 'data' structure */
    gmx_qsort(data, rotg->nat, sizeof(sort_along_vec_t), projection_compare);
<<<<<<< HEAD
    
    /* Copy back the sorted values */
    for (i=0; i<rotg->nat; i++)
    {
        copy_rvec(data[i].x    , erg->xc[i]           );
=======

    /* Copy back the sorted values */
    for (i = 0; i < rotg->nat; i++)
    {
        copy_rvec(data[i].x, erg->xc[i]           );
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        copy_rvec(data[i].x_ref, erg->xc_ref_sorted[i]);
        erg->mc_sorted[i]  = data[i].m;
        erg->xc_sortind[i] = data[i].ind;
    }
}


/* For each slab, get the first and the last index of the sorted atom
 * indices */
static void get_firstlast_atom_per_slab(t_rotgrp *rotg)
{
<<<<<<< HEAD
    int i,islab,n;
    real beta;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /* Find the first atom that needs to enter the calculation for each slab */
    n = erg->slab_first;  /* slab */
    i = 0; /* start with the first atom */
=======
    int             i, islab, n;
    real            beta;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;

    /* Find the first atom that needs to enter the calculation for each slab */
    n = erg->slab_first; /* slab */
    i = 0;               /* start with the first atom */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    do
    {
        /* Find the first atom that significantly contributes to this slab */
        do /* move forward in position until a large enough beta is found */
        {
            beta = calc_beta(erg->xc[i], rotg, n);
            i++;
<<<<<<< HEAD
        } while ((beta < -erg->max_beta) && (i < rotg->nat));
        i--;
        islab = n - erg->slab_first;  /* slab index */
        erg->firstatom[islab] = i;
        /* Proceed to the next slab */
        n++;
    } while (n <= erg->slab_last);
    
    /* Find the last atom for each slab */
     n = erg->slab_last; /* start with last slab */
     i = rotg->nat-1;  /* start with the last atom */
     do
     {
         do /* move backward in position until a large enough beta is found */
         {
             beta = calc_beta(erg->xc[i], rotg, n);
             i--;
         } while ((beta > erg->max_beta) && (i > -1));
         i++;
         islab = n - erg->slab_first;  /* slab index */
         erg->lastatom[islab] = i;
         /* Proceed to the next slab */
         n--;
     } while (n >= erg->slab_first);
}


/* Determine the very first and very last slab that needs to be considered 
 * For the first slab that needs to be considered, we have to find the smallest
 * n that obeys:
 * 
 * x_first * v - n*Delta_x <= beta_max
 * 
 * slab index n, slab distance Delta_x, rotation vector v. For the last slab we 
 * have to find the largest n that obeys
 * 
 * x_last * v - n*Delta_x >= -beta_max
 *  
 */
static gmx_inline int get_first_slab(
        t_rotgrp *rotg,     /* The rotation group (inputrec data) */
        real     max_beta,  /* The max_beta value, instead of min_gaussian */
        rvec     firstatom) /* First atom after sorting along the rotation vector v */
{
    /* Find the first slab for the first atom */   
=======
        }
        while ((beta < -erg->max_beta) && (i < rotg->nat));
        i--;
        islab                 = n - erg->slab_first; /* slab index */
        erg->firstatom[islab] = i;
        /* Proceed to the next slab */
        n++;
    }
    while (n <= erg->slab_last);

    /* Find the last atom for each slab */
    n = erg->slab_last; /* start with last slab */
    i = rotg->nat-1;    /* start with the last atom */
    do
    {
        do  /* move backward in position until a large enough beta is found */
        {
            beta = calc_beta(erg->xc[i], rotg, n);
            i--;
        }
        while ((beta > erg->max_beta) && (i > -1));
        i++;
        islab                = n - erg->slab_first; /* slab index */
        erg->lastatom[islab] = i;
        /* Proceed to the next slab */
        n--;
    }
    while (n >= erg->slab_first);
}


/* Determine the very first and very last slab that needs to be considered
 * For the first slab that needs to be considered, we have to find the smallest
 * n that obeys:
 *
 * x_first * v - n*Delta_x <= beta_max
 *
 * slab index n, slab distance Delta_x, rotation vector v. For the last slab we
 * have to find the largest n that obeys
 *
 * x_last * v - n*Delta_x >= -beta_max
 *
 */
static gmx_inline int get_first_slab(
        t_rotgrp *rotg,      /* The rotation group (inputrec data) */
        real      max_beta,  /* The max_beta value, instead of min_gaussian */
        rvec      firstatom) /* First atom after sorting along the rotation vector v */
{
    /* Find the first slab for the first atom */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return ceil((iprod(firstatom, rotg->vec) - max_beta)/rotg->slab_dist);
}


static gmx_inline int get_last_slab(
        t_rotgrp *rotg,     /* The rotation group (inputrec data) */
<<<<<<< HEAD
        real     max_beta,  /* The max_beta value, instead of min_gaussian */
        rvec     lastatom)  /* Last atom along v */
{
    /* Find the last slab for the last atom */
    return floor((iprod(lastatom, rotg->vec) + max_beta)/rotg->slab_dist);    
=======
        real      max_beta, /* The max_beta value, instead of min_gaussian */
        rvec      lastatom) /* Last atom along v */
{
    /* Find the last slab for the last atom */
    return floor((iprod(lastatom, rotg->vec) + max_beta)/rotg->slab_dist);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


static void get_firstlast_slab_check(
<<<<<<< HEAD
        t_rotgrp        *rotg,     /* The rotation group (inputrec data) */
        t_gmx_enfrotgrp *erg,      /* The rotation group (data only accessible in this file) */
        rvec            firstatom, /* First atom after sorting along the rotation vector v */
        rvec            lastatom,  /* Last atom along v */
        int             g)         /* The rotation group number */
=======
        t_rotgrp        *rotg,      /* The rotation group (inputrec data) */
        t_gmx_enfrotgrp *erg,       /* The rotation group (data only accessible in this file) */
        rvec             firstatom, /* First atom after sorting along the rotation vector v */
        rvec             lastatom,  /* Last atom along v */
        int              g)         /* The rotation group number */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    erg->slab_first = get_first_slab(rotg, erg->max_beta, firstatom);
    erg->slab_last  = get_last_slab(rotg, erg->max_beta, lastatom);

    /* Check whether we have reference data to compare against */
    if (erg->slab_first < erg->slab_first_ref)
<<<<<<< HEAD
        gmx_fatal(FARGS, "%s No reference data for first slab (n=%d), unable to proceed.",
                  RotStr, erg->slab_first);
    
    /* Check whether we have reference data to compare against */
    if (erg->slab_last > erg->slab_last_ref)
        gmx_fatal(FARGS, "%s No reference data for last slab (n=%d), unable to proceed.",
                  RotStr, erg->slab_last);
=======
    {
        gmx_fatal(FARGS, "%s No reference data for first slab (n=%d), unable to proceed.",
                  RotStr, erg->slab_first);
    }

    /* Check whether we have reference data to compare against */
    if (erg->slab_last > erg->slab_last_ref)
    {
        gmx_fatal(FARGS, "%s No reference data for last slab (n=%d), unable to proceed.",
                  RotStr, erg->slab_last);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* Enforced rotation with a flexible axis */
static void do_flexible(
<<<<<<< HEAD
        gmx_bool  bMaster,
        gmx_enfrot_t enfrot,    /* Other rotation data                        */
        t_rotgrp  *rotg,        /* The rotation group                         */
        int       g,            /* Group number                               */
        rvec      x[],          /* The local positions                        */
        matrix    box,
        double    t,            /* Time in picoseconds                        */
        gmx_large_int_t step,   /* The time step                              */
        gmx_bool  bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool  bOutstepSlab) /* Output per-slab data                       */
{
    int          l,nslabs;
    real         sigma;       /* The Gaussian width sigma */
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;
    
    /* Sort the collective coordinates erg->xc along the rotation vector. This is
     * an optimization for the inner loop. */
    sort_collective_coordinates(rotg, enfrot->data);
    
    /* Determine the first relevant slab for the first atom and the last
     * relevant slab for the last atom */
    get_firstlast_slab_check(rotg, erg, erg->xc[0], erg->xc[rotg->nat-1], g);
    
=======
        gmx_bool        bMaster,
        gmx_enfrot_t    enfrot,       /* Other rotation data                        */
        t_rotgrp       *rotg,         /* The rotation group                         */
        int             g,            /* Group number                               */
        rvec            x[],          /* The local positions                        */
        matrix          box,
        double          t,            /* Time in picoseconds                        */
        gmx_large_int_t step,         /* The time step                              */
        gmx_bool        bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool        bOutstepSlab) /* Output per-slab data                       */
{
    int             l, nslabs;
    real            sigma;    /* The Gaussian width sigma */
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */


    erg = rotg->enfrotgrp;

    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;

    /* Sort the collective coordinates erg->xc along the rotation vector. This is
     * an optimization for the inner loop. */
    sort_collective_coordinates(rotg, enfrot->data);

    /* Determine the first relevant slab for the first atom and the last
     * relevant slab for the last atom */
    get_firstlast_slab_check(rotg, erg, erg->xc[0], erg->xc[rotg->nat-1], g);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Determine for each slab depending on the min_gaussian cutoff criterium,
     * a first and a last atom index inbetween stuff needs to be calculated */
    get_firstlast_atom_per_slab(rotg);

    /* Determine the gaussian-weighted center of positions for all slabs */
<<<<<<< HEAD
    get_slab_centers(rotg,erg->xc,erg->mc_sorted,g,t,enfrot->out_slabs,bOutstepSlab,FALSE);
        
    /* Clear the torque per slab from last time step: */
    nslabs = erg->slab_last - erg->slab_first + 1;
    for (l=0; l<nslabs; l++)
        erg->slab_torque_v[l] = 0.0;
    
    /* Call the rotational forces kernel */
    if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEXT)
        erg->V = do_flex_lowlevel(rotg, sigma, x, bOutstepRot, bOutstepSlab, box);
    else if (rotg->eType == erotgFLEX2 || rotg->eType == erotgFLEX2T)
        erg->V = do_flex2_lowlevel(rotg, sigma, x, bOutstepRot, bOutstepSlab, box);
    else
        gmx_fatal(FARGS, "Unknown flexible rotation type");
    
    /* Determine angle by RMSD fit to the reference - Let's hope this */
    /* only happens once in a while, since this is not parallelized! */
    if ( bMaster && (erotgFitPOT != rotg->eFittype) )
=======
    get_slab_centers(rotg, erg->xc, erg->mc_sorted, g, t, enfrot->out_slabs, bOutstepSlab, FALSE);

    /* Clear the torque per slab from last time step: */
    nslabs = erg->slab_last - erg->slab_first + 1;
    for (l = 0; l < nslabs; l++)
    {
        erg->slab_torque_v[l] = 0.0;
    }

    /* Call the rotational forces kernel */
    if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEXT)
    {
        erg->V = do_flex_lowlevel(rotg, sigma, x, bOutstepRot, bOutstepSlab, box);
    }
    else if (rotg->eType == erotgFLEX2 || rotg->eType == erotgFLEX2T)
    {
        erg->V = do_flex2_lowlevel(rotg, sigma, x, bOutstepRot, bOutstepSlab, box);
    }
    else
    {
        gmx_fatal(FARGS, "Unknown flexible rotation type");
    }

    /* Determine angle by RMSD fit to the reference - Let's hope this */
    /* only happens once in a while, since this is not parallelized! */
    if (bMaster && (erotgFitPOT != rotg->eFittype) )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (bOutstepRot)
        {
            /* Fit angle of the whole rotation group */
            erg->angle_v = flex_fit_angle(rotg);
        }
        if (bOutstepSlab)
        {
            /* Fit angle of each slab */
            flex_fit_angle_perslab(g, rotg, t, erg->degangle, enfrot->out_angles);
        }
    }

    /* Lump together the torques from all slabs: */
    erg->torque_v = 0.0;
<<<<<<< HEAD
    for (l=0; l<nslabs; l++)
         erg->torque_v += erg->slab_torque_v[l];
=======
    for (l = 0; l < nslabs; l++)
    {
        erg->torque_v += erg->slab_torque_v[l];
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* Calculate the angle between reference and actual rotation group atom,
 * both projected into a plane perpendicular to the rotation vector: */
static void angle(t_rotgrp *rotg,
<<<<<<< HEAD
        rvec x_act,
        rvec x_ref,
        real *alpha,
        real *weight)  /* atoms near the rotation axis should count less than atoms far away */
{
    rvec xp, xrp;  /* current and reference positions projected on a plane perpendicular to pg->vec */
=======
                  rvec      x_act,
                  rvec      x_ref,
                  real     *alpha,
                  real     *weight) /* atoms near the rotation axis should count less than atoms far away */
{
    rvec xp, xrp;                   /* current and reference positions projected on a plane perpendicular to pg->vec */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rvec dum;


    /* Project x_ref and x into a plane through the origin perpendicular to rot_vec: */
    /* Project x_ref: xrp = x_ref - (vec * x_ref) * vec */
    svmul(iprod(rotg->vec, x_ref), rotg->vec, dum);
    rvec_sub(x_ref, dum, xrp);
    /* Project x_act: */
    svmul(iprod(rotg->vec, x_act), rotg->vec, dum);
    rvec_sub(x_act, dum, xp);

    /* Retrieve information about which vector precedes. gmx_angle always
     * returns a positive angle. */
    cprod(xp, xrp, dum); /* if reference precedes, this is pointing into the same direction as vec */

    if (iprod(rotg->vec, dum) >= 0)
<<<<<<< HEAD
        *alpha = -gmx_angle(xrp, xp);
    else
        *alpha = +gmx_angle(xrp, xp);
    
=======
    {
        *alpha = -gmx_angle(xrp, xp);
    }
    else
    {
        *alpha = +gmx_angle(xrp, xp);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Also return the weight */
    *weight = norm(xp);
}


<<<<<<< HEAD
/* Project first vector onto a plane perpendicular to the second vector 
=======
/* Project first vector onto a plane perpendicular to the second vector
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * dr = dr - (dr.v)v
 * Note that v must be of unit length.
 */
static gmx_inline void project_onto_plane(rvec dr, const rvec v)
{
    rvec tmp;
<<<<<<< HEAD
    
    
    svmul(iprod(dr,v),v,tmp);  /* tmp = (dr.v)v */
    rvec_dec(dr, tmp);         /* dr = dr - (dr.v)v */
=======


    svmul(iprod(dr, v), v, tmp); /* tmp = (dr.v)v */
    rvec_dec(dr, tmp);           /* dr = dr - (dr.v)v */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* Fixed rotation: The rotation reference group rotates around the v axis. */
/* The atoms of the actual rotation group are attached with imaginary  */
/* springs to the reference atoms.                                     */
static void do_fixed(
<<<<<<< HEAD
        t_rotgrp  *rotg,        /* The rotation group                         */
        rvec      x[],          /* The positions                              */
        matrix    box,          /* The simulation box                         */
        double    t,            /* Time in picoseconds                        */
        gmx_large_int_t step,   /* The time step                              */
        gmx_bool  bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool  bOutstepSlab) /* Output per-slab data                       */
{
    int       ifit,j,jj,m;
    rvec      dr;
    rvec      tmp_f;           /* Force */
    real      alpha;           /* a single angle between an actual and a reference position */
    real      weight;          /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec      xi_xc;           /* xi - xc */
    gmx_bool  bCalcPotFit;
    rvec      fit_xr_loc;
=======
        t_rotgrp       *rotg,         /* The rotation group                         */
        rvec            x[],          /* The positions                              */
        matrix          box,          /* The simulation box                         */
        double          t,            /* Time in picoseconds                        */
        gmx_large_int_t step,         /* The time step                              */
        gmx_bool        bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool        bOutstepSlab) /* Output per-slab data                       */
{
    int             ifit, j, jj, m;
    rvec            dr;
    rvec            tmp_f;     /* Force */
    real            alpha;     /* a single angle between an actual and a reference position */
    real            weight;    /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec            xi_xc;     /* xi - xc */
    gmx_bool        bCalcPotFit;
    rvec            fit_xr_loc;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* for mass weighting: */
    real      wi;              /* Mass-weighting of the positions */
    real      N_M;             /* N/M */
    real      k_wi;            /* k times wi */

    gmx_bool  bProject;

<<<<<<< HEAD
    
    erg=rotg->enfrotgrp;
    bProject = (rotg->eType==erotgPM) || (rotg->eType==erotgPMPF);
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT==rotg->eFittype);
=======

    erg         = rotg->enfrotgrp;
    bProject    = (rotg->eType == erotgPM) || (rotg->eType == erotgPMPF);
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT == rotg->eFittype);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    N_M = rotg->nat * erg->invmass;

    /* Each process calculates the forces on its local atoms */
<<<<<<< HEAD
    for (j=0; j<erg->nat_loc; j++)
=======
    for (j = 0; j < erg->nat_loc; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Calculate (x_i-x_c) resp. (x_i-u) */
        rvec_sub(erg->x_loc_pbc[j], erg->xc_center, xi_xc);

        /* Calculate Omega*(y_i-y_c)-(x_i-x_c) */
        rvec_sub(erg->xr_loc[j], xi_xc, dr);
<<<<<<< HEAD
        
        if (bProject)
            project_onto_plane(dr, rotg->vec);
            
=======

        if (bProject)
        {
            project_onto_plane(dr, rotg->vec);
        }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Mass-weighting */
        wi = N_M*erg->m_loc[j];

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        k_wi = rotg->k*wi;
<<<<<<< HEAD
        for (m=0; m<DIM; m++)
=======
        for (m = 0; m < DIM; m++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            tmp_f[m]             = k_wi*dr[m];
            erg->f_rot_loc[j][m] = tmp_f[m];
            erg->V              += 0.5*k_wi*sqr(dr[m]);
        }
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (ifit = 0; ifit < rotg->PotAngle_nstep; ifit++)
            {
                /* Index of this rotation group atom with respect to the whole rotation group */
                jj = erg->xc_ref_ind[j];

                /* Rotate with the alternative angle. Like rotate_local_reference(),
                 * just for a single local atom */
                mvmul(erg->PotAngleFit->rotmat[ifit], rotg->x_ref[jj], fit_xr_loc); /* fit_xr_loc = Omega*(y_i-y_c) */

                /* Calculate Omega*(y_i-y_c)-(x_i-x_c) */
                rvec_sub(fit_xr_loc, xi_xc, dr);

                if (bProject)
<<<<<<< HEAD
                    project_onto_plane(dr, rotg->vec);
=======
                {
                    project_onto_plane(dr, rotg->vec);
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5*k_wi*norm2(dr);
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(rotg->vec, tmp_f, erg->x_loc_pbc[j], erg->xc_center);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Calculate the angle between reference and actual rotation group atom. */
            angle(rotg, xi_xc, erg->xr_loc[j], &alpha, &weight);  /* angle in rad, weighted */
            erg->angle_v  += alpha * weight;
            erg->weight_v += weight;
        }
        /* If you want enforced rotation to contribute to the virial,
         * activate the following lines:
            if (MASTER(cr))
            {
               Add the rotation contribution to the virial
              for(j=0; j<DIM; j++)
                for(m=0;m<DIM;m++)
                  vir[j][m] += 0.5*f[ii][j]*dr[m];
            }
         */

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
}


/* Calculate the radial motion potential and forces */
static void do_radial_motion(
<<<<<<< HEAD
        t_rotgrp  *rotg,        /* The rotation group                         */
        rvec      x[],          /* The positions                              */
        matrix    box,          /* The simulation box                         */
        double    t,            /* Time in picoseconds                        */
        gmx_large_int_t step,   /* The time step                              */
        gmx_bool  bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool  bOutstepSlab) /* Output per-slab data                       */
{
    int       j,jj,ifit;
    rvec      tmp_f;           /* Force */
    real      alpha;           /* a single angle between an actual and a reference position */
    real      weight;          /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec      xj_u;            /* xj - u */
    rvec      tmpvec,fit_tmpvec;
    real      fac,fac2,sum=0.0;
    rvec      pj;
    gmx_bool  bCalcPotFit;
=======
        t_rotgrp       *rotg,         /* The rotation group                         */
        rvec            x[],          /* The positions                              */
        matrix          box,          /* The simulation box                         */
        double          t,            /* Time in picoseconds                        */
        gmx_large_int_t step,         /* The time step                              */
        gmx_bool        bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool        bOutstepSlab) /* Output per-slab data                       */
{
    int             j, jj, ifit;
    rvec            tmp_f;     /* Force */
    real            alpha;     /* a single angle between an actual and a reference position */
    real            weight;    /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec            xj_u;      /* xj - u */
    rvec            tmpvec, fit_tmpvec;
    real            fac, fac2, sum = 0.0;
    rvec            pj;
    gmx_bool        bCalcPotFit;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* For mass weighting: */
    real      wj;              /* Mass-weighting of the positions */
    real      N_M;             /* N/M */


<<<<<<< HEAD
    erg=rotg->enfrotgrp;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT==rotg->eFittype);
=======
    erg         = rotg->enfrotgrp;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT == rotg->eFittype);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    N_M = rotg->nat * erg->invmass;

    /* Each process calculates the forces on its local atoms */
<<<<<<< HEAD
    for (j=0; j<erg->nat_loc; j++)
=======
    for (j = 0; j < erg->nat_loc; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Calculate (xj-u) */
        rvec_sub(erg->x_loc_pbc[j], erg->xc_center, xj_u);  /* xj_u = xj-u */

        /* Calculate Omega.(yj0-u) */
        cprod(rotg->vec, erg->xr_loc[j], tmpvec);  /* tmpvec = v x Omega.(yj0-u) */

<<<<<<< HEAD
                              /*         v x Omega.(yj0-u)     */
        unitv(tmpvec, pj);    /*  pj = ---------------------   */
                              /*       | v x Omega.(yj0-u) |   */

        fac = iprod(pj, xj_u);  /* fac = pj.(xj-u) */
=======
        /*                       *         v x Omega.(yj0-u)     */
        unitv(tmpvec, pj);      /*  pj = ---------------------   */
                                /*       | v x Omega.(yj0-u) |   */

        fac  = iprod(pj, xj_u); /* fac = pj.(xj-u) */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fac2 = fac*fac;

        /* Mass-weighting */
        wj = N_M*erg->m_loc[j];

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        svmul(-rotg->k*wj*fac, pj, tmp_f);
        copy_rvec(tmp_f, erg->f_rot_loc[j]);
        sum += wj*fac2;

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (ifit = 0; ifit < rotg->PotAngle_nstep; ifit++)
            {
                /* Index of this rotation group atom with respect to the whole rotation group */
                jj = erg->xc_ref_ind[j];

                /* Rotate with the alternative angle. Like rotate_local_reference(),
                 * just for a single local atom */
                mvmul(erg->PotAngleFit->rotmat[ifit], rotg->x_ref[jj], fit_tmpvec); /* fit_tmpvec = Omega*(yj0-u) */

                /* Calculate Omega.(yj0-u) */
<<<<<<< HEAD
                cprod(rotg->vec, fit_tmpvec, tmpvec);  /* tmpvec = v x Omega.(yj0-u) */
                                      /*         v x Omega.(yj0-u)     */
                unitv(tmpvec, pj);    /*  pj = ---------------------   */
                                      /*       | v x Omega.(yj0-u) |   */

                fac = iprod(pj, xj_u);  /* fac = pj.(xj-u) */
=======
                cprod(rotg->vec, fit_tmpvec, tmpvec); /* tmpvec = v x Omega.(yj0-u) */
                /*                                     *         v x Omega.(yj0-u)     */
                unitv(tmpvec, pj);                    /*  pj = ---------------------   */
                                                      /*       | v x Omega.(yj0-u) |   */

                fac  = iprod(pj, xj_u);               /* fac = pj.(xj-u) */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                fac2 = fac*fac;

                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5*rotg->k*wj*fac2;
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(rotg->vec, tmp_f, erg->x_loc_pbc[j], erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(rotg, xj_u, erg->xr_loc[j], &alpha, &weight);  /* angle in rad, weighted */
            erg->angle_v  += alpha * weight;
            erg->weight_v += weight;
        }

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
    erg->V = 0.5*rotg->k*sum;
}


/* Calculate the radial motion pivot-free potential and forces */
static void do_radial_motion_pf(
<<<<<<< HEAD
        t_rotgrp  *rotg,        /* The rotation group                         */
        rvec      x[],          /* The positions                              */
        matrix    box,          /* The simulation box                         */
        double    t,            /* Time in picoseconds                        */
        gmx_large_int_t step,   /* The time step                              */
        gmx_bool  bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool  bOutstepSlab) /* Output per-slab data                       */
{
    int       i,ii,iigrp,ifit,j;
    rvec      xj;              /* Current position */
    rvec      xj_xc;           /* xj  - xc  */
    rvec      yj0_yc0;         /* yj0 - yc0 */
    rvec      tmp_f;           /* Force */
    real      alpha;           /* a single angle between an actual and a reference position */
    real      weight;          /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec      tmpvec, tmpvec2;
    rvec      innersumvec;     /* Precalculation of the inner sum */
    rvec      innersumveckM;
    real      fac,fac2,V=0.0;
    rvec      qi,qj;
    gmx_bool  bCalcPotFit;

    /* For mass weighting: */
    real      mj,wi,wj;        /* Mass-weighting of the positions */
    real      N_M;             /* N/M */


    erg=rotg->enfrotgrp;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT==rotg->eFittype);
=======
        t_rotgrp       *rotg,         /* The rotation group                         */
        rvec            x[],          /* The positions                              */
        matrix          box,          /* The simulation box                         */
        double          t,            /* Time in picoseconds                        */
        gmx_large_int_t step,         /* The time step                              */
        gmx_bool        bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool        bOutstepSlab) /* Output per-slab data                       */
{
    int             i, ii, iigrp, ifit, j;
    rvec            xj;          /* Current position */
    rvec            xj_xc;       /* xj  - xc  */
    rvec            yj0_yc0;     /* yj0 - yc0 */
    rvec            tmp_f;       /* Force */
    real            alpha;       /* a single angle between an actual and a reference position */
    real            weight;      /* single weight for a single angle */
    gmx_enfrotgrp_t erg;         /* Pointer to enforced rotation group data */
    rvec            tmpvec, tmpvec2;
    rvec            innersumvec; /* Precalculation of the inner sum */
    rvec            innersumveckM;
    real            fac, fac2, V = 0.0;
    rvec            qi, qj;
    gmx_bool        bCalcPotFit;

    /* For mass weighting: */
    real      mj, wi, wj;      /* Mass-weighting of the positions */
    real      N_M;             /* N/M */


    erg         = rotg->enfrotgrp;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT == rotg->eFittype);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    N_M = rotg->nat * erg->invmass;

    /* Get the current center of the rotation group: */
    get_center(erg->xc, erg->mc, rotg->nat, erg->xc_center);

    /* Precalculate Sum_i [ wi qi.(xi-xc) qi ] which is needed for every single j */
    clear_rvec(innersumvec);
<<<<<<< HEAD
    for (i=0; i < rotg->nat; i++)
=======
    for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Mass-weighting */
        wi = N_M*erg->mc[i];

        /* Calculate qi. Note that xc_ref_center has already been subtracted from
         * x_ref in init_rot_group.*/
<<<<<<< HEAD
        mvmul(erg->rotmat, rotg->x_ref[i], tmpvec);  /* tmpvec  = Omega.(yi0-yc0) */

        cprod(rotg->vec, tmpvec, tmpvec2);          /* tmpvec2 = v x Omega.(yi0-yc0) */

                              /*         v x Omega.(yi0-yc0)     */
        unitv(tmpvec2, qi);   /*  qi = -----------------------   */
                              /*       | v x Omega.(yi0-yc0) |   */

        rvec_sub(erg->xc[i], erg->xc_center, tmpvec);  /* tmpvec = xi-xc */
=======
        mvmul(erg->rotmat, rotg->x_ref[i], tmpvec); /* tmpvec  = Omega.(yi0-yc0) */

        cprod(rotg->vec, tmpvec, tmpvec2);          /* tmpvec2 = v x Omega.(yi0-yc0) */

        /*                                             *         v x Omega.(yi0-yc0)     */
        unitv(tmpvec2, qi);                           /*  qi = -----------------------   */
                                                      /*       | v x Omega.(yi0-yc0) |   */

        rvec_sub(erg->xc[i], erg->xc_center, tmpvec); /* tmpvec = xi-xc */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        svmul(wi*iprod(qi, tmpvec), qi, tmpvec2);

        rvec_inc(innersumvec, tmpvec2);
    }
    svmul(rotg->k*erg->invmass, innersumvec, innersumveckM);

    /* Each process calculates the forces on its local atoms */
<<<<<<< HEAD
    for (j=0; j<erg->nat_loc; j++)
=======
    for (j = 0; j < erg->nat_loc; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Local index of a rotation group atom  */
        ii = erg->ind_loc[j];
        /* Position of this atom in the collective array */
        iigrp = erg->xc_ref_ind[j];
        /* Mass-weighting */
        mj = erg->mc[iigrp];  /* need the unsorted mass here */
        wj = N_M*mj;

        /* Current position of this atom: x[ii][XX/YY/ZZ] */
        copy_rvec(x[ii], xj);

        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

        /* The (unrotated) reference position is yj0. yc0 has already
         * been subtracted in init_rot_group */
        copy_rvec(rotg->x_ref[iigrp], yj0_yc0);   /* yj0_yc0 = yj0 - yc0      */

        /* Calculate Omega.(yj0-yc0) */
<<<<<<< HEAD
        mvmul(erg->rotmat, yj0_yc0, tmpvec2);     /* tmpvec2 = Omega.(yj0 - yc0)  */

        cprod(rotg->vec, tmpvec2, tmpvec);  /* tmpvec = v x Omega.(yj0-yc0) */

                              /*         v x Omega.(yj0-yc0)     */
=======
        mvmul(erg->rotmat, yj0_yc0, tmpvec2); /* tmpvec2 = Omega.(yj0 - yc0)  */

        cprod(rotg->vec, tmpvec2, tmpvec);    /* tmpvec = v x Omega.(yj0-yc0) */

        /*                     *         v x Omega.(yj0-yc0)     */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        unitv(tmpvec, qj);    /*  qj = -----------------------   */
                              /*       | v x Omega.(yj0-yc0) |   */

        /* Calculate (xj-xc) */
<<<<<<< HEAD
        rvec_sub(xj, erg->xc_center, xj_xc);  /* xj_xc = xj-xc */

        fac = iprod(qj, xj_xc);  /* fac = qj.(xj-xc) */
=======
        rvec_sub(xj, erg->xc_center, xj_xc); /* xj_xc = xj-xc */

        fac  = iprod(qj, xj_xc);             /* fac = qj.(xj-xc) */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fac2 = fac*fac;

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        svmul(-rotg->k*wj*fac, qj, tmp_f); /* part 1 of force */
        svmul(mj, innersumveckM, tmpvec);  /* part 2 of force */
        rvec_inc(tmp_f, tmpvec);
        copy_rvec(tmp_f, erg->f_rot_loc[j]);
        V += wj*fac2;

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (ifit = 0; ifit < rotg->PotAngle_nstep; ifit++)
            {
                /* Rotate with the alternative angle. Like rotate_local_reference(),
                 * just for a single local atom */
                mvmul(erg->PotAngleFit->rotmat[ifit], yj0_yc0, tmpvec2); /* tmpvec2 = Omega*(yj0-yc0) */

                /* Calculate Omega.(yj0-u) */
<<<<<<< HEAD
                cprod(rotg->vec, tmpvec2, tmpvec);  /* tmpvec = v x Omega.(yj0-yc0) */
                                      /*         v x Omega.(yj0-yc0)     */
                unitv(tmpvec, qj);    /*  qj = -----------------------   */
                                      /*       | v x Omega.(yj0-yc0) |   */

                fac = iprod(qj, xj_xc);  /* fac = qj.(xj-xc) */
=======
                cprod(rotg->vec, tmpvec2, tmpvec); /* tmpvec = v x Omega.(yj0-yc0) */
                /*                                  *         v x Omega.(yj0-yc0)     */
                unitv(tmpvec, qj);                 /*  qj = -----------------------   */
                                                   /*       | v x Omega.(yj0-yc0) |   */

                fac  = iprod(qj, xj_xc);           /* fac = qj.(xj-xc) */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                fac2 = fac*fac;

                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5*rotg->k*wj*fac2;
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(rotg->vec, tmp_f, xj, erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(rotg, xj_xc, yj0_yc0, &alpha, &weight);  /* angle in rad, weighted */
            erg->angle_v  += alpha * weight;
            erg->weight_v += weight;
        }

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
    erg->V = 0.5*rotg->k*V;
}


/* Precalculate the inner sum for the radial motion 2 forces */
static void radial_motion2_precalc_inner_sum(t_rotgrp  *rotg, rvec innersumvec)
{
<<<<<<< HEAD
    int       i;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec      xi_xc;           /* xj - xc */
    rvec      tmpvec,tmpvec2;
    real      fac,fac2;
    rvec      ri,si;
    real      siri;
    rvec      v_xi_xc;          /* v x (xj - u) */
    real      psii,psiistar;
    real      wi;              /* Mass-weighting of the positions */
    real      N_M;             /* N/M */
    rvec      sumvec;

    erg=rotg->enfrotgrp;
=======
    int             i;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec            xi_xc;     /* xj - xc */
    rvec            tmpvec, tmpvec2;
    real            fac, fac2;
    rvec            ri, si;
    real            siri;
    rvec            v_xi_xc;   /* v x (xj - u) */
    real            psii, psiistar;
    real            wi;        /* Mass-weighting of the positions */
    real            N_M;       /* N/M */
    rvec            sumvec;

    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    N_M = rotg->nat * erg->invmass;

    /* Loop over the collective set of positions */
    clear_rvec(sumvec);
<<<<<<< HEAD
    for (i=0; i<rotg->nat; i++)
=======
    for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Mass-weighting */
        wi = N_M*erg->mc[i];

        rvec_sub(erg->xc[i], erg->xc_center, xi_xc); /* xi_xc = xi-xc         */

        /* Calculate ri. Note that xc_ref_center has already been subtracted from
         * x_ref in init_rot_group.*/
        mvmul(erg->rotmat, rotg->x_ref[i], ri);      /* ri  = Omega.(yi0-yc0) */

        cprod(rotg->vec, xi_xc, v_xi_xc);            /* v_xi_xc = v x (xi-u)  */

        fac = norm2(v_xi_xc);
<<<<<<< HEAD
                                          /*                      1           */
=======
        /*                                 *                      1           */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        psiistar = 1.0/(fac + rotg->eps); /* psiistar = --------------------- */
                                          /*            |v x (xi-xc)|^2 + eps */

        psii = gmx_invsqrt(fac);          /*                 1                */
                                          /*  psii    = -------------         */
                                          /*            |v x (xi-xc)|         */

<<<<<<< HEAD
        svmul(psii, v_xi_xc, si);          /*  si = psii * (v x (xi-xc) )     */

        fac = iprod(v_xi_xc, ri);                   /* fac = (v x (xi-xc)).ri */
=======
        svmul(psii, v_xi_xc, si);         /*  si = psii * (v x (xi-xc) )     */

        fac  = iprod(v_xi_xc, ri);        /* fac = (v x (xi-xc)).ri */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fac2 = fac*fac;

        siri = iprod(si, ri);                       /* siri = si.ri           */

        svmul(psiistar/psii, ri, tmpvec);
        svmul(psiistar*psiistar/(psii*psii*psii) * siri, si, tmpvec2);
        rvec_dec(tmpvec, tmpvec2);
        cprod(tmpvec, rotg->vec, tmpvec2);

        svmul(wi*siri, tmpvec2, tmpvec);

        rvec_inc(sumvec, tmpvec);
    }
    svmul(rotg->k*erg->invmass, sumvec, innersumvec);
}


/* Calculate the radial motion 2 potential and forces */
static void do_radial_motion2(
<<<<<<< HEAD
        t_rotgrp  *rotg,        /* The rotation group                         */
        rvec      x[],          /* The positions                              */
        matrix    box,          /* The simulation box                         */
        double    t,            /* Time in picoseconds                        */
        gmx_large_int_t step,   /* The time step                              */
        gmx_bool  bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool  bOutstepSlab) /* Output per-slab data                       */
{
    int       ii,iigrp,ifit,j;
    rvec      xj;              /* Position */
    real      alpha;           /* a single angle between an actual and a reference position */
    real      weight;          /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec      xj_u;            /* xj - u */
    rvec      yj0_yc0;         /* yj0 -yc0 */
    rvec      tmpvec,tmpvec2;
    real      fac,fit_fac,fac2,Vpart=0.0;
    rvec      rj,fit_rj,sj;
    real      sjrj;
    rvec      v_xj_u;          /* v x (xj - u) */
    real      psij,psijstar;
    real      mj,wj;           /* For mass-weighting of the positions */
    real      N_M;             /* N/M */
    gmx_bool  bPF;
    rvec      innersumvec;
    gmx_bool  bCalcPotFit;


    erg=rotg->enfrotgrp;

    bPF = rotg->eType==erotgRM2PF;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT==rotg->eFittype);
=======
        t_rotgrp       *rotg,         /* The rotation group                         */
        rvec            x[],          /* The positions                              */
        matrix          box,          /* The simulation box                         */
        double          t,            /* Time in picoseconds                        */
        gmx_large_int_t step,         /* The time step                              */
        gmx_bool        bOutstepRot,  /* Output to main rotation output file        */
        gmx_bool        bOutstepSlab) /* Output per-slab data                       */
{
    int             ii, iigrp, ifit, j;
    rvec            xj;        /* Position */
    real            alpha;     /* a single angle between an actual and a reference position */
    real            weight;    /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec            xj_u;      /* xj - u */
    rvec            yj0_yc0;   /* yj0 -yc0 */
    rvec            tmpvec, tmpvec2;
    real            fac, fit_fac, fac2, Vpart = 0.0;
    rvec            rj, fit_rj, sj;
    real            sjrj;
    rvec            v_xj_u;    /* v x (xj - u) */
    real            psij, psijstar;
    real            mj, wj;    /* For mass-weighting of the positions */
    real            N_M;       /* N/M */
    gmx_bool        bPF;
    rvec            innersumvec;
    gmx_bool        bCalcPotFit;


    erg = rotg->enfrotgrp;

    bPF         = rotg->eType == erotgRM2PF;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (erotgFitPOT == rotg->eFittype);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    clear_rvec(yj0_yc0); /* Make the compiler happy */

    clear_rvec(innersumvec);
    if (bPF)
    {
        /* For the pivot-free variant we have to use the current center of
         * mass of the rotation group instead of the pivot u */
        get_center(erg->xc, erg->mc, rotg->nat, erg->xc_center);

        /* Also, we precalculate the second term of the forces that is identical
         * (up to the weight factor mj) for all forces */
<<<<<<< HEAD
        radial_motion2_precalc_inner_sum(rotg,innersumvec);
=======
        radial_motion2_precalc_inner_sum(rotg, innersumvec);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    N_M = rotg->nat * erg->invmass;

    /* Each process calculates the forces on its local atoms */
<<<<<<< HEAD
    for (j=0; j<erg->nat_loc; j++)
=======
    for (j = 0; j < erg->nat_loc; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (bPF)
        {
            /* Local index of a rotation group atom  */
            ii = erg->ind_loc[j];
            /* Position of this atom in the collective array */
            iigrp = erg->xc_ref_ind[j];
            /* Mass-weighting */
            mj = erg->mc[iigrp];

            /* Current position of this atom: x[ii] */
            copy_rvec(x[ii], xj);

            /* Shift this atom such that it is near its reference */
            shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

            /* The (unrotated) reference position is yj0. yc0 has already
             * been subtracted in init_rot_group */
            copy_rvec(rotg->x_ref[iigrp], yj0_yc0);   /* yj0_yc0 = yj0 - yc0  */

            /* Calculate Omega.(yj0-yc0) */
            mvmul(erg->rotmat, yj0_yc0, rj);         /* rj = Omega.(yj0-yc0)  */
        }
        else
        {
            mj = erg->m_loc[j];
            copy_rvec(erg->x_loc_pbc[j], xj);
            copy_rvec(erg->xr_loc[j], rj);           /* rj = Omega.(yj0-u)    */
        }
        /* Mass-weighting */
        wj = N_M*mj;

        /* Calculate (xj-u) resp. (xj-xc) */
        rvec_sub(xj, erg->xc_center, xj_u);          /* xj_u = xj-u           */

        cprod(rotg->vec, xj_u, v_xj_u);              /* v_xj_u = v x (xj-u)   */

        fac = norm2(v_xj_u);
<<<<<<< HEAD
                                          /*                      1           */
=======
        /*                                 *                      1           */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        psijstar = 1.0/(fac + rotg->eps); /*  psistar = --------------------  */
                                          /*            |v x (xj-u)|^2 + eps  */

        psij = gmx_invsqrt(fac);          /*                 1                */
                                          /*  psij    = ------------          */
                                          /*            |v x (xj-u)|          */

        svmul(psij, v_xj_u, sj);          /*  sj = psij * (v x (xj-u) )       */

<<<<<<< HEAD
        fac = iprod(v_xj_u, rj);                     /* fac = (v x (xj-u)).rj */
=======
        fac  = iprod(v_xj_u, rj);         /* fac = (v x (xj-u)).rj */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fac2 = fac*fac;

        sjrj = iprod(sj, rj);                        /* sjrj = sj.rj          */

        svmul(psijstar/psij, rj, tmpvec);
        svmul(psijstar*psijstar/(psij*psij*psij) * sjrj, sj, tmpvec2);
        rvec_dec(tmpvec, tmpvec2);
        cprod(tmpvec, rotg->vec, tmpvec2);

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        svmul(-rotg->k*wj*sjrj, tmpvec2, tmpvec);
        svmul(mj, innersumvec, tmpvec2);  /* This is != 0 only for the pivot-free variant */

        rvec_add(tmpvec2, tmpvec, erg->f_rot_loc[j]);
        Vpart += wj*psijstar*fac2;

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (ifit = 0; ifit < rotg->PotAngle_nstep; ifit++)
            {
                if (bPF)
                {
                    mvmul(erg->PotAngleFit->rotmat[ifit], yj0_yc0, fit_rj); /* fit_rj = Omega.(yj0-yc0) */
                }
                else
                {
                    /* Position of this atom in the collective array */
                    iigrp = erg->xc_ref_ind[j];
                    /* Rotate with the alternative angle. Like rotate_local_reference(),
                     * just for a single local atom */
                    mvmul(erg->PotAngleFit->rotmat[ifit], rotg->x_ref[iigrp], fit_rj); /* fit_rj = Omega*(yj0-u) */
                }
<<<<<<< HEAD
                fit_fac = iprod(v_xj_u, fit_rj); /* fac = (v x (xj-u)).fit_rj */
=======
                fit_fac = iprod(v_xj_u, fit_rj);                                       /* fac = (v x (xj-u)).fit_rj */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5*rotg->k*wj*psijstar*fit_fac*fit_fac;
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(rotg->vec, erg->f_rot_loc[j], xj, erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(rotg, xj_u, rj, &alpha, &weight);  /* angle in rad, weighted */
            erg->angle_v  += alpha * weight;
            erg->weight_v += weight;
        }

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
    erg->V = 0.5*rotg->k*Vpart;
}


<<<<<<< HEAD
/* Determine the smallest and largest position vector (with respect to the 
 * rotation vector) for the reference group */
static void get_firstlast_atom_ref(
        t_rotgrp  *rotg, 
        int       *firstindex, 
        int       *lastindex)
{
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    int i;
    real xcproj;               /* The projection of a reference position on the 
                                  rotation vector */
    real minproj, maxproj;     /* Smallest and largest projection on v */
    

    
    erg=rotg->enfrotgrp;
=======
/* Determine the smallest and largest position vector (with respect to the
 * rotation vector) for the reference group */
static void get_firstlast_atom_ref(
        t_rotgrp  *rotg,
        int       *firstindex,
        int       *lastindex)
{
    gmx_enfrotgrp_t erg;              /* Pointer to enforced rotation group data */
    int             i;
    real            xcproj;           /* The projection of a reference position on the
                                         rotation vector */
    real            minproj, maxproj; /* Smallest and largest projection on v */



    erg = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Start with some value */
    minproj = iprod(rotg->x_ref[0], rotg->vec);
    maxproj = minproj;
<<<<<<< HEAD
    
    /* This is just to ensure that it still works if all the atoms of the 
     * reference structure are situated in a plane perpendicular to the rotation 
     * vector */
    *firstindex = 0;
    *lastindex  = rotg->nat-1;
    
    /* Loop over all atoms of the reference group, 
     * project them on the rotation vector to find the extremes */
    for (i=0; i<rotg->nat; i++)
=======

    /* This is just to ensure that it still works if all the atoms of the
     * reference structure are situated in a plane perpendicular to the rotation
     * vector */
    *firstindex = 0;
    *lastindex  = rotg->nat-1;

    /* Loop over all atoms of the reference group,
     * project them on the rotation vector to find the extremes */
    for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        xcproj = iprod(rotg->x_ref[i], rotg->vec);
        if (xcproj < minproj)
        {
<<<<<<< HEAD
            minproj = xcproj;
=======
            minproj     = xcproj;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            *firstindex = i;
        }
        if (xcproj > maxproj)
        {
<<<<<<< HEAD
            maxproj = xcproj;
=======
            maxproj    = xcproj;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            *lastindex = i;
        }
    }
}


/* Allocate memory for the slabs */
static void allocate_slabs(
<<<<<<< HEAD
        t_rotgrp  *rotg, 
        FILE      *fplog, 
        int       g, 
        gmx_bool  bVerbose)
{
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    int i, nslabs;
    
    
    erg=rotg->enfrotgrp;
    
    /* More slabs than are defined for the reference are never needed */
    nslabs = erg->slab_last_ref - erg->slab_first_ref + 1;
    
    /* Remember how many we allocated */
    erg->nslabs_alloc = nslabs;

    if ( (NULL != fplog) && bVerbose )
        fprintf(fplog, "%s allocating memory to store data for %d slabs (rotation group %d).\n",
                RotStr, nslabs,g);
    snew(erg->slab_center     , nslabs);
    snew(erg->slab_center_ref , nslabs);
    snew(erg->slab_weights    , nslabs);
    snew(erg->slab_torque_v   , nslabs);
    snew(erg->slab_data       , nslabs);
    snew(erg->gn_atom         , nslabs);
    snew(erg->gn_slabind      , nslabs);
    snew(erg->slab_innersumvec, nslabs);
    for (i=0; i<nslabs; i++)
    {
        snew(erg->slab_data[i].x     , rotg->nat);
        snew(erg->slab_data[i].ref   , rotg->nat);
        snew(erg->slab_data[i].weight, rotg->nat);
    }
    snew(erg->xc_ref_sorted, rotg->nat);
    snew(erg->xc_sortind   , rotg->nat);
    snew(erg->firstatom    , nslabs);
    snew(erg->lastatom     , nslabs);
}


/* From the extreme coordinates of the reference group, determine the first 
=======
        t_rotgrp  *rotg,
        FILE      *fplog,
        int        g,
        gmx_bool   bVerbose)
{
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    int             i, nslabs;


    erg = rotg->enfrotgrp;

    /* More slabs than are defined for the reference are never needed */
    nslabs = erg->slab_last_ref - erg->slab_first_ref + 1;

    /* Remember how many we allocated */
    erg->nslabs_alloc = nslabs;

    if ( (NULL != fplog) && bVerbose)
    {
        fprintf(fplog, "%s allocating memory to store data for %d slabs (rotation group %d).\n",
                RotStr, nslabs, g);
    }
    snew(erg->slab_center, nslabs);
    snew(erg->slab_center_ref, nslabs);
    snew(erg->slab_weights, nslabs);
    snew(erg->slab_torque_v, nslabs);
    snew(erg->slab_data, nslabs);
    snew(erg->gn_atom, nslabs);
    snew(erg->gn_slabind, nslabs);
    snew(erg->slab_innersumvec, nslabs);
    for (i = 0; i < nslabs; i++)
    {
        snew(erg->slab_data[i].x, rotg->nat);
        snew(erg->slab_data[i].ref, rotg->nat);
        snew(erg->slab_data[i].weight, rotg->nat);
    }
    snew(erg->xc_ref_sorted, rotg->nat);
    snew(erg->xc_sortind, rotg->nat);
    snew(erg->firstatom, nslabs);
    snew(erg->lastatom, nslabs);
}


/* From the extreme coordinates of the reference group, determine the first
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * and last slab of the reference. We can never have more slabs in the real
 * simulation than calculated here for the reference.
 */
static void get_firstlast_slab_ref(t_rotgrp *rotg, real mc[], int ref_firstindex, int ref_lastindex)
{
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
<<<<<<< HEAD
    int first,last,firststart;
    rvec dummy;

    
    erg=rotg->enfrotgrp;
    first = get_first_slab(rotg, erg->max_beta, rotg->x_ref[ref_firstindex]);
    last  = get_last_slab( rotg, erg->max_beta, rotg->x_ref[ref_lastindex ]);
=======
    int             first, last, firststart;
    rvec            dummy;


    erg        = rotg->enfrotgrp;
    first      = get_first_slab(rotg, erg->max_beta, rotg->x_ref[ref_firstindex]);
    last       = get_last_slab( rotg, erg->max_beta, rotg->x_ref[ref_lastindex ]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    firststart = first;

    while (get_slab_weight(first, rotg, rotg->x_ref, mc, &dummy) > WEIGHT_MIN)
    {
        first--;
    }
    erg->slab_first_ref = first+1;
    while (get_slab_weight(last, rotg, rotg->x_ref, mc, &dummy) > WEIGHT_MIN)
    {
        last++;
    }
    erg->slab_last_ref  = last-1;
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    erg->slab_buffer = firststart - erg->slab_first_ref;
}


<<<<<<< HEAD

static void init_rot_group(FILE *fplog,t_commrec *cr,int g,t_rotgrp *rotg,
        rvec *x,gmx_mtop_t *mtop,gmx_bool bVerbose,FILE *out_slabs, gmx_bool bOutputCenters)
{
    int i,ii;
    rvec        coord,*xdum;
    gmx_bool    bFlex,bColl;
    t_atom      *atom;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    int         ref_firstindex, ref_lastindex;
    real        mass,totalmass;
    real        start=0.0;
    
=======
/* Special version of copy_rvec:
 * During the copy procedure of xcurr to b, the correct PBC image is chosen
 * such that the copied vector ends up near its reference position xref */
static inline void copy_correct_pbc_image(
        const rvec  xcurr,  /* copy vector xcurr ...                */
        rvec        b,      /* ... to b ...                         */
        const rvec  xref,   /* choosing the PBC image such that b ends up near xref */
        matrix      box,
        int         npbcdim)
{
    rvec  dx;
    int   d, m;
    ivec  shift;


    /* Shortest PBC distance between the atom and its reference */
    rvec_sub(xcurr, xref, dx);

    /* Determine the shift for this atom */
    clear_ivec(shift);
    for (m = npbcdim-1; m >= 0; m--)
    {
        while (dx[m] < -0.5*box[m][m])
        {
            for (d = 0; d < DIM; d++)
            {
                dx[d] += box[m][d];
            }
            shift[m]++;
        }
        while (dx[m] >= 0.5*box[m][m])
        {
            for (d = 0; d < DIM; d++)
            {
                dx[d] -= box[m][d];
            }
            shift[m]--;
        }
    }

    /* Apply the shift to the position */
    copy_rvec(xcurr, b);
    shift_single_coord(box, b, shift);
}


static void init_rot_group(FILE *fplog, t_commrec *cr, int g, t_rotgrp *rotg,
                           rvec *x, gmx_mtop_t *mtop, gmx_bool bVerbose, FILE *out_slabs, matrix box,
                           gmx_bool bOutputCenters)
{
    int                   i, ii;
    rvec                  coord, *xdum;
    gmx_bool              bFlex, bColl;
    t_atom               *atom;
    gmx_enfrotgrp_t       erg; /* Pointer to enforced rotation group data */
    int                   ref_firstindex, ref_lastindex;
    gmx_mtop_atomlookup_t alook = NULL;
    real                  mass, totalmass;
    real                  start = 0.0;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Do we have a flexible axis? */
    bFlex = ISFLEX(rotg);
    /* Do we use a global set of coordinates? */
    bColl = ISCOLL(rotg);

<<<<<<< HEAD
    erg=rotg->enfrotgrp;
    
    /* Allocate space for collective coordinates if needed */
    if (bColl)
    {
        snew(erg->xc        , rotg->nat);
        snew(erg->xc_shifts , rotg->nat);
=======
    erg = rotg->enfrotgrp;

    /* Allocate space for collective coordinates if needed */
    if (bColl)
    {
        snew(erg->xc, rotg->nat);
        snew(erg->xc_shifts, rotg->nat);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        snew(erg->xc_eshifts, rotg->nat);

        /* Save the original (whole) set of positions such that later the
         * molecule can always be made whole again */
<<<<<<< HEAD
        snew(erg->xc_old    , rotg->nat);        
        if (MASTER(cr))
        {
            for (i=0; i<rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], erg->xc_old[i]);
=======
        snew(erg->xc_old, rotg->nat);
        if (MASTER(cr))
        {
            for (i = 0; i < rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_correct_pbc_image(x[ii], erg->xc_old[i], rotg->x_ref[i], box, 3);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
#ifdef GMX_MPI
        if (PAR(cr))
<<<<<<< HEAD
            gmx_bcast(rotg->nat*sizeof(erg->xc_old[0]),erg->xc_old, cr);
=======
        {
            gmx_bcast(rotg->nat*sizeof(erg->xc_old[0]), erg->xc_old, cr);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

        if (rotg->eFittype == erotgFitNORM)
        {
            snew(erg->xc_ref_length, rotg->nat); /* in case fit type NORM is chosen */
<<<<<<< HEAD
            snew(erg->xc_norm      , rotg->nat);
=======
            snew(erg->xc_norm, rotg->nat);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
    else
    {
<<<<<<< HEAD
        snew(erg->xr_loc   , rotg->nat);
        snew(erg->x_loc_pbc, rotg->nat);
    }
    
    snew(erg->f_rot_loc , rotg->nat);
    snew(erg->xc_ref_ind, rotg->nat);
    
=======
        snew(erg->xr_loc, rotg->nat);
        snew(erg->x_loc_pbc, rotg->nat);
    }

    snew(erg->f_rot_loc, rotg->nat);
    snew(erg->xc_ref_ind, rotg->nat);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Make space for the calculation of the potential at other angles (used
     * for fitting only) */
    if (erotgFitPOT == rotg->eFittype)
    {
        snew(erg->PotAngleFit, 1);
        snew(erg->PotAngleFit->degangle, rotg->PotAngle_nstep);
<<<<<<< HEAD
        snew(erg->PotAngleFit->V       , rotg->PotAngle_nstep);
        snew(erg->PotAngleFit->rotmat  , rotg->PotAngle_nstep);
=======
        snew(erg->PotAngleFit->V, rotg->PotAngle_nstep);
        snew(erg->PotAngleFit->rotmat, rotg->PotAngle_nstep);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* Get the set of angles around the reference angle */
        start = -0.5 * (rotg->PotAngle_nstep - 1)*rotg->PotAngle_step;
        for (i = 0; i < rotg->PotAngle_nstep; i++)
<<<<<<< HEAD
            erg->PotAngleFit->degangle[i] = start + i*rotg->PotAngle_step;
=======
        {
            erg->PotAngleFit->degangle[i] = start + i*rotg->PotAngle_step;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        erg->PotAngleFit = NULL;
    }

    /* xc_ref_ind needs to be set to identity in the serial case */
    if (!PAR(cr))
<<<<<<< HEAD
        for (i=0; i<rotg->nat; i++)
            erg->xc_ref_ind[i] = i;

    /* Copy the masses so that the center can be determined. For all types of
     * enforced rotation, we store the masses in the erg->mc array. */
    snew(erg->mc, rotg->nat);
    if (bFlex)
        snew(erg->mc_sorted, rotg->nat);
    if (!bColl)
        snew(erg->m_loc, rotg->nat);
    totalmass=0.0;
    for (i=0; i<rotg->nat; i++)
    {
        if (rotg->bMassW)
        {
            gmx_mtop_atomnr_to_atom(mtop,rotg->ind[i],&atom);
            mass=atom->m;
        }
        else
        {
            mass=1.0;
=======
    {
        for (i = 0; i < rotg->nat; i++)
        {
            erg->xc_ref_ind[i] = i;
        }
    }

    /* Copy the masses so that the center can be determined. For all types of
     * enforced rotation, we store the masses in the erg->mc array. */
    if (rotg->bMassW)
    {
        alook = gmx_mtop_atomlookup_init(mtop);
    }
    snew(erg->mc, rotg->nat);
    if (bFlex)
    {
        snew(erg->mc_sorted, rotg->nat);
    }
    if (!bColl)
    {
        snew(erg->m_loc, rotg->nat);
    }
    totalmass = 0.0;
    for (i = 0; i < rotg->nat; i++)
    {
        if (rotg->bMassW)
        {
            gmx_mtop_atomnr_to_atom(alook, rotg->ind[i], &atom);
            mass = atom->m;
        }
        else
        {
            mass = 1.0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        erg->mc[i] = mass;
        totalmass += mass;
    }
    erg->invmass = 1.0/totalmass;
<<<<<<< HEAD
    
    /* Set xc_ref_center for any rotation potential */
    if ((rotg->eType==erotgISO) || (rotg->eType==erotgPM) || (rotg->eType==erotgRM) || (rotg->eType==erotgRM2))
=======

    if (rotg->bMassW)
    {
        gmx_mtop_atomlookup_destroy(alook);
    }

    /* Set xc_ref_center for any rotation potential */
    if ((rotg->eType == erotgISO) || (rotg->eType == erotgPM) || (rotg->eType == erotgRM) || (rotg->eType == erotgRM2))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Set the pivot point for the fixed, stationary-axis potentials. This
         * won't change during the simulation */
        copy_rvec(rotg->pivot, erg->xc_ref_center);
        copy_rvec(rotg->pivot, erg->xc_center    );
    }
    else
    {
        /* Center of the reference positions */
        get_center(rotg->x_ref, erg->mc, rotg->nat, erg->xc_ref_center);

        /* Center of the actual positions */
        if (MASTER(cr))
        {
            snew(xdum, rotg->nat);
<<<<<<< HEAD
            for (i=0; i<rotg->nat; i++)
=======
            for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], xdum[i]);
            }
            get_center(xdum, erg->mc, rotg->nat, erg->xc_center);
            sfree(xdum);
        }
#ifdef GMX_MPI
        if (PAR(cr))
<<<<<<< HEAD
            gmx_bcast(sizeof(erg->xc_center), erg->xc_center, cr);
=======
        {
            gmx_bcast(sizeof(erg->xc_center), erg->xc_center, cr);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    }

    if ( (rotg->eType != erotgFLEX) && (rotg->eType != erotgFLEX2) )
    {
        /* Put the reference positions into origin: */
<<<<<<< HEAD
        for (i=0; i<rotg->nat; i++)
            rvec_dec(rotg->x_ref[i], erg->xc_ref_center);
=======
        for (i = 0; i < rotg->nat; i++)
        {
            rvec_dec(rotg->x_ref[i], erg->xc_ref_center);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    /* Enforced rotation with flexible axis */
    if (bFlex)
    {
        /* Calculate maximum beta value from minimum gaussian (performance opt.) */
        erg->max_beta = calc_beta_max(rotg->min_gaussian, rotg->slab_dist);

        /* Determine the smallest and largest coordinate with respect to the rotation vector */
        get_firstlast_atom_ref(rotg, &ref_firstindex, &ref_lastindex);
<<<<<<< HEAD
        
        /* From the extreme coordinates of the reference group, determine the first 
         * and last slab of the reference. */
        get_firstlast_slab_ref(rotg, erg->mc, ref_firstindex, ref_lastindex);
                
=======

        /* From the extreme coordinates of the reference group, determine the first
         * and last slab of the reference. */
        get_firstlast_slab_ref(rotg, erg->mc, ref_firstindex, ref_lastindex);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Allocate memory for the slabs */
        allocate_slabs(rotg, fplog, g, bVerbose);

        /* Flexible rotation: determine the reference centers for the rest of the simulation */
        erg->slab_first = erg->slab_first_ref;
<<<<<<< HEAD
        erg->slab_last = erg->slab_last_ref;
        get_slab_centers(rotg,rotg->x_ref,erg->mc,g,-1,out_slabs,bOutputCenters,TRUE);
=======
        erg->slab_last  = erg->slab_last_ref;
        get_slab_centers(rotg, rotg->x_ref, erg->mc, g, -1, out_slabs, bOutputCenters, TRUE);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* Length of each x_rotref vector from center (needed if fit routine NORM is chosen): */
        if (rotg->eFittype == erotgFitNORM)
        {
<<<<<<< HEAD
            for (i=0; i<rotg->nat; i++)
=======
            for (i = 0; i < rotg->nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                rvec_sub(rotg->x_ref[i], erg->xc_ref_center, coord);
                erg->xc_ref_length[i] = norm(coord);
            }
        }
    }
}


<<<<<<< HEAD
extern void dd_make_local_rotation_groups(gmx_domdec_t *dd,t_rot *rot)
{
    gmx_ga2la_t ga2la;
    int g;
    t_rotgrp *rotg;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    
    ga2la = dd->ga2la;

    for(g=0; g<rot->ngrp; g++)
=======
extern void dd_make_local_rotation_groups(gmx_domdec_t *dd, t_rot *rot)
{
    gmx_ga2la_t     ga2la;
    int             g;
    t_rotgrp       *rotg;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */

    ga2la = dd->ga2la;

    for (g = 0; g < rot->ngrp; g++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        rotg = &rot->grp[g];
        erg  = rotg->enfrotgrp;


<<<<<<< HEAD
        dd_make_local_group_indices(ga2la,rotg->nat,rotg->ind,
                &erg->nat_loc,&erg->ind_loc,&erg->nalloc_loc,erg->xc_ref_ind);
=======
        dd_make_local_group_indices(ga2la, rotg->nat, rotg->ind,
                                    &erg->nat_loc, &erg->ind_loc, &erg->nalloc_loc, erg->xc_ref_ind);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}


/* Calculate the size of the MPI buffer needed in reduce_output() */
static int calc_mpi_bufsize(t_rot *rot)
{
<<<<<<< HEAD
    int g;
    int count_group, count_total;
    t_rotgrp *rotg;
=======
    int             g;
    int             count_group, count_total;
    t_rotgrp       *rotg;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */


    count_total = 0;
<<<<<<< HEAD
    for (g=0; g<rot->ngrp; g++)
=======
    for (g = 0; g < rot->ngrp; g++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        rotg = &rot->grp[g];
        erg  = rotg->enfrotgrp;

        /* Count the items that are transferred for this group: */
        count_group = 4; /* V, torque, angle, weight */

        /* Add the maximum number of slabs for flexible groups */
        if (ISFLEX(rotg))
<<<<<<< HEAD
            count_group += erg->slab_last_ref - erg->slab_first_ref + 1;

        /* Add space for the potentials at different angles: */
        if (erotgFitPOT == rotg->eFittype)
            count_group += rotg->PotAngle_nstep;
=======
        {
            count_group += erg->slab_last_ref - erg->slab_first_ref + 1;
        }

        /* Add space for the potentials at different angles: */
        if (erotgFitPOT == rotg->eFittype)
        {
            count_group += rotg->PotAngle_nstep;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* Add to the total number: */
        count_total += count_group;
    }

    return count_total;
}


<<<<<<< HEAD
extern void init_rot(FILE *fplog,t_inputrec *ir,int nfile,const t_filenm fnm[],
        t_commrec *cr, rvec *x, matrix box, gmx_mtop_t *mtop, const output_env_t oenv,
        gmx_bool bVerbose, unsigned long Flags)
{
    t_rot    *rot;
    t_rotgrp *rotg;
    int      g;
    int      nat_max=0;     /* Size of biggest rotation group */
    gmx_enfrot_t er;        /* Pointer to the enforced rotation buffer variables */    
    gmx_enfrotgrp_t erg;    /* Pointer to enforced rotation group data */
    rvec     *x_pbc=NULL;   /* Space for the pbc-correct atom positions */


    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
        gmx_fatal(FARGS, "Enforced rotation is only implemented for domain decomposition!");

    if ( MASTER(cr) && bVerbose)
        fprintf(stdout, "%s Initializing ...\n", RotStr);

    rot = ir->rot;
    snew(rot->enfrot, 1);
    er = rot->enfrot;
=======
extern void init_rot(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
                     t_commrec *cr, rvec *x, matrix box, gmx_mtop_t *mtop, const output_env_t oenv,
                     gmx_bool bVerbose, unsigned long Flags)
{
    t_rot          *rot;
    t_rotgrp       *rotg;
    int             g;
    int             nat_max = 0;  /* Size of biggest rotation group */
    gmx_enfrot_t    er;           /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg;          /* Pointer to enforced rotation group data */
    rvec           *x_pbc = NULL; /* Space for the pbc-correct atom positions */


    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
    {
        gmx_fatal(FARGS, "Enforced rotation is only implemented for domain decomposition!");
    }

    if (MASTER(cr) && bVerbose)
    {
        fprintf(stdout, "%s Initializing ...\n", RotStr);
    }

    rot = ir->rot;
    snew(rot->enfrot, 1);
    er        = rot->enfrot;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    er->Flags = Flags;

    /* When appending, skip first output to avoid duplicate entries in the data files */
    if (er->Flags & MD_APPENDFILES)
<<<<<<< HEAD
        er->bOut = FALSE;
    else
        er->bOut = TRUE;

    if ( MASTER(cr) && er->bOut )
        please_cite(fplog, "Kutzner2011");
=======
    {
        er->bOut = FALSE;
    }
    else
    {
        er->bOut = TRUE;
    }

    if (MASTER(cr) && er->bOut)
    {
        please_cite(fplog, "Kutzner2011");
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Output every step for reruns */
    if (er->Flags & MD_RERUN)
    {
        if (NULL != fplog)
<<<<<<< HEAD
            fprintf(fplog, "%s rerun - will write rotation output every available step.\n", RotStr);
=======
        {
            fprintf(fplog, "%s rerun - will write rotation output every available step.\n", RotStr);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        rot->nstrout = 1;
        rot->nstsout = 1;
    }

    er->out_slabs = NULL;
<<<<<<< HEAD
    if ( MASTER(cr) && HaveFlexibleGroups(rot) )
        er->out_slabs = open_slab_out(opt2fn("-rs",nfile,fnm), rot, oenv);
=======
    if (MASTER(cr) && HaveFlexibleGroups(rot) )
    {
        er->out_slabs = open_slab_out(opt2fn("-rs", nfile, fnm), rot, oenv);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (MASTER(cr))
    {
        /* Remove pbc, make molecule whole.
         * When ir->bContinuation=TRUE this has already been done, but ok. */
<<<<<<< HEAD
        snew(x_pbc,mtop->natoms);
        m_rveccopy(mtop->natoms,x,x_pbc);
        do_pbc_first_mtop(NULL,ir->ePBC,box,mtop,x_pbc);
    }

    for (g=0; g<rot->ngrp; g++)
=======
        snew(x_pbc, mtop->natoms);
        m_rveccopy(mtop->natoms, x, x_pbc);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        /* All molecules will be whole now, but not necessarily in the home box.
         * Additionally, if a rotation group consists of more than one molecule
         * (e.g. two strands of DNA), each one of them can end up in a different
         * periodic box. This is taken care of in init_rot_group.  */
    }

    for (g = 0; g < rot->ngrp; g++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        rotg = &rot->grp[g];

        if (NULL != fplog)
<<<<<<< HEAD
            fprintf(fplog,"%s group %d type '%s'\n", RotStr, g, erotg_names[rotg->eType]);
        
=======
        {
            fprintf(fplog, "%s group %d type '%s'\n", RotStr, g, erotg_names[rotg->eType]);
        }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (rotg->nat > 0)
        {
            /* Allocate space for the rotation group's data: */
            snew(rotg->enfrotgrp, 1);
            erg  = rotg->enfrotgrp;

<<<<<<< HEAD
            nat_max=max(nat_max, rotg->nat);
            
=======
            nat_max = max(nat_max, rotg->nat);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (PAR(cr))
            {
                erg->nat_loc    = 0;
                erg->nalloc_loc = 0;
                erg->ind_loc    = NULL;
            }
            else
            {
                erg->nat_loc = rotg->nat;
                erg->ind_loc = rotg->ind;
            }
<<<<<<< HEAD
            init_rot_group(fplog,cr,g,rotg,x_pbc,mtop,bVerbose,er->out_slabs,
=======
            init_rot_group(fplog, cr, g, rotg, x_pbc, mtop, bVerbose, er->out_slabs, box,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                           !(er->Flags & MD_APPENDFILES) ); /* Do not output the reference centers
                                                             * again if we are appending */
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Allocate space for enforced rotation buffer variables */
    er->bufsize = nat_max;
    snew(er->data, nat_max);
    snew(er->xbuf, nat_max);
    snew(er->mbuf, nat_max);

    /* Buffers for MPI reducing torques, angles, weights (for each group), and V */
    if (PAR(cr))
    {
        er->mpi_bufsize = calc_mpi_bufsize(rot) + 100; /* larger to catch errors */
<<<<<<< HEAD
        snew(er->mpi_inbuf , er->mpi_bufsize);
=======
        snew(er->mpi_inbuf, er->mpi_bufsize);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        snew(er->mpi_outbuf, er->mpi_bufsize);
    }
    else
    {
        er->mpi_bufsize = 0;
<<<<<<< HEAD
        er->mpi_inbuf = NULL;
        er->mpi_outbuf = NULL;
=======
        er->mpi_inbuf   = NULL;
        er->mpi_outbuf  = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    /* Only do I/O on the MASTER */
    er->out_angles  = NULL;
    er->out_rot     = NULL;
    er->out_torque  = NULL;
    if (MASTER(cr))
    {
<<<<<<< HEAD
        er->out_rot = open_rot_out(opt2fn("-ro",nfile,fnm), rot, oenv);

        if (rot->nstsout > 0)
        {
            if ( HaveFlexibleGroups(rot) || HavePotFitGroups(rot) )
                er->out_angles  = open_angles_out(opt2fn("-ra",nfile,fnm), rot, oenv);
            if ( HaveFlexibleGroups(rot) )
                er->out_torque  = open_torque_out(opt2fn("-rt",nfile,fnm), rot, oenv);
=======
        er->out_rot = open_rot_out(opt2fn("-ro", nfile, fnm), rot, oenv);

        if (rot->nstsout > 0)
        {
            if (HaveFlexibleGroups(rot) || HavePotFitGroups(rot) )
            {
                er->out_angles  = open_angles_out(opt2fn("-ra", nfile, fnm), rot, oenv);
            }
            if (HaveFlexibleGroups(rot) )
            {
                er->out_torque  = open_torque_out(opt2fn("-rt", nfile, fnm), rot, oenv);
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }

        sfree(x_pbc);
    }
}


<<<<<<< HEAD
extern void finish_rot(FILE *fplog,t_rot *rot)
{
    gmx_enfrot_t er;        /* Pointer to the enforced rotation buffer variables */    

    
    er=rot->enfrot;
    if (er->out_rot)
        gmx_fio_fclose(er->out_rot);
    if (er->out_slabs)
        gmx_fio_fclose(er->out_slabs);
    if (er->out_angles)
        gmx_fio_fclose(er->out_angles);
    if (er->out_torque)
        gmx_fio_fclose(er->out_torque);
=======
extern void finish_rot(FILE *fplog, t_rot *rot)
{
    gmx_enfrot_t er;        /* Pointer to the enforced rotation buffer variables */


    er = rot->enfrot;
    if (er->out_rot)
    {
        gmx_fio_fclose(er->out_rot);
    }
    if (er->out_slabs)
    {
        gmx_fio_fclose(er->out_slabs);
    }
    if (er->out_angles)
    {
        gmx_fio_fclose(er->out_angles);
    }
    if (er->out_torque)
    {
        gmx_fio_fclose(er->out_torque);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* Rotate the local reference positions and store them in
 * erg->xr_loc[0...(nat_loc-1)]
 *
 * Note that we already subtracted u or y_c from the reference positions
 * in init_rot_group().
 */
static void rotate_local_reference(t_rotgrp *rotg)
{
    gmx_enfrotgrp_t erg;
<<<<<<< HEAD
    int i,ii;

    
    erg=rotg->enfrotgrp;
    
    for (i=0; i<erg->nat_loc; i++)
=======
    int             i, ii;


    erg = rotg->enfrotgrp;

    for (i = 0; i < erg->nat_loc; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Index of this rotation group atom with respect to the whole rotation group */
        ii = erg->xc_ref_ind[i];
        /* Rotate */
        mvmul(erg->rotmat, rotg->x_ref[ii], erg->xr_loc[i]);
    }
}


/* Select the PBC representation for each local x position and store that
 * for later usage. We assume the right PBC image of an x is the one nearest to
 * its rotated reference */
static void choose_pbc_image(rvec x[], t_rotgrp *rotg, matrix box, int npbcdim)
{
<<<<<<< HEAD
    int d,i,ii,m;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec xref,xcurr,dx;
    ivec shift;


    erg=rotg->enfrotgrp;
    
    for (i=0; i<erg->nat_loc; i++)
    {
        clear_ivec(shift);
        
=======
    int             i, ii;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec            xref;


    erg = rotg->enfrotgrp;

    for (i = 0; i < erg->nat_loc; i++)
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Index of a rotation group atom  */
        ii = erg->ind_loc[i];

        /* Get the reference position. The pivot was already
         * subtracted in init_rot_group() from the reference positions. Also,
         * the reference positions have already been rotated in
         * rotate_local_reference() */
        copy_rvec(erg->xr_loc[i], xref);
<<<<<<< HEAD
        
        /* Subtract the (old) center from the current positions
         * (just to determine the shifts!) */
        rvec_sub(x[ii], erg->xc_center, xcurr);
        
        /* Shortest PBC distance between the atom and its reference */
        rvec_sub(xcurr, xref, dx);
        
        /* Determine the shift for this atom */
        for(m=npbcdim-1; m>=0; m--)
        {
            while (dx[m] < -0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] += box[m][d];
                shift[m]++;
            }
            while (dx[m] >= 0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] -= box[m][d];
                shift[m]--;
            }
        }
        
        /* Apply the shift to the current atom */
        copy_rvec(x[ii], erg->x_loc_pbc[i]);
        shift_single_coord(box, erg->x_loc_pbc[i], shift); 
=======

        copy_correct_pbc_image(x[ii], erg->x_loc_pbc[i], xref, box, npbcdim);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}


extern void do_rotation(
<<<<<<< HEAD
        t_commrec *cr,
        t_inputrec *ir,
        matrix box,
        rvec x[],
        real t,
        gmx_large_int_t step,
        gmx_wallcycle_t wcycle,
        gmx_bool bNS)
{
    int      g,i,ii;
    t_rot    *rot;
    t_rotgrp *rotg;
    gmx_bool outstep_slab, outstep_rot;
    gmx_bool bFlex,bColl;
    gmx_enfrot_t er;     /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg; /* Pointer to enforced rotation group data           */
    rvec     transvec;
    t_gmx_potfit *fit=NULL;     /* For fit type 'potential' determine the fit
=======
        t_commrec      *cr,
        t_inputrec     *ir,
        matrix          box,
        rvec            x[],
        real            t,
        gmx_large_int_t step,
        gmx_wallcycle_t wcycle,
        gmx_bool        bNS)
{
    int             g, i, ii;
    t_rot          *rot;
    t_rotgrp       *rotg;
    gmx_bool        outstep_slab, outstep_rot;
    gmx_bool        bFlex, bColl;
    gmx_enfrot_t    er;         /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg;        /* Pointer to enforced rotation group data           */
    rvec            transvec;
    t_gmx_potfit   *fit = NULL; /* For fit type 'potential' determine the fit
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                   angle via the potential minimum            */

    /* Enforced rotation cycle counting: */
    gmx_cycles_t cycles_comp;   /* Cycles for the enf. rotation computation
                                   only, does not count communication. This
                                   counter is used for load-balancing         */

#ifdef TAKETIME
    double t0;
#endif
<<<<<<< HEAD
    
    rot=ir->rot;
    er=rot->enfrot;
    
=======

    rot = ir->rot;
    er  = rot->enfrot;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* When to output in main rotation output file */
    outstep_rot  = do_per_step(step, rot->nstrout) && er->bOut;
    /* When to output per-slab data */
    outstep_slab = do_per_step(step, rot->nstsout) && er->bOut;

    /* Output time into rotation output file */
    if (outstep_rot && MASTER(cr))
<<<<<<< HEAD
        fprintf(er->out_rot, "%12.3e",t);

    /**************************************************************************/
    /* First do ALL the communication! */
    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg=rotg->enfrotgrp;
=======
    {
        fprintf(er->out_rot, "%12.3e", t);
    }

    /**************************************************************************/
    /* First do ALL the communication! */
    for (g = 0; g < rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg  = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* Do we have a flexible axis? */
        bFlex = ISFLEX(rotg);
        /* Do we use a collective (global) set of coordinates? */
        bColl = ISCOLL(rotg);

        /* Calculate the rotation matrix for this angle: */
        erg->degangle = rotg->rate * t;
<<<<<<< HEAD
        calc_rotmat(rotg->vec,erg->degangle,erg->rotmat);
=======
        calc_rotmat(rotg->vec, erg->degangle, erg->rotmat);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        if (bColl)
        {
            /* Transfer the rotation group's positions such that every node has
             * all of them. Every node contributes its local positions x and stores
             * it in the collective erg->xc array. */
<<<<<<< HEAD
            communicate_group_positions(cr,erg->xc, erg->xc_shifts, erg->xc_eshifts, bNS,
                    x, rotg->nat, erg->nat_loc, erg->ind_loc, erg->xc_ref_ind, erg->xc_old, box);
=======
            communicate_group_positions(cr, erg->xc, erg->xc_shifts, erg->xc_eshifts, bNS,
                                        x, rotg->nat, erg->nat_loc, erg->ind_loc, erg->xc_ref_ind, erg->xc_old, box);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        else
        {
            /* Fill the local masses array;
             * this array changes in DD/neighborsearching steps */
            if (bNS)
            {
<<<<<<< HEAD
                for (i=0; i<erg->nat_loc; i++)
                {
                    /* Index of local atom w.r.t. the collective rotation group */
                    ii = erg->xc_ref_ind[i];
=======
                for (i = 0; i < erg->nat_loc; i++)
                {
                    /* Index of local atom w.r.t. the collective rotation group */
                    ii            = erg->xc_ref_ind[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    erg->m_loc[i] = erg->mc[ii];
                }
            }

            /* Calculate Omega*(y_i-y_c) for the local positions */
            rotate_local_reference(rotg);

            /* Choose the nearest PBC images of the group atoms with respect
             * to the rotated reference positions */
            choose_pbc_image(x, rotg, box, 3);

            /* Get the center of the rotation group */
<<<<<<< HEAD
            if ( (rotg->eType==erotgISOPF) || (rotg->eType==erotgPMPF) )
                get_center_comm(cr, erg->x_loc_pbc, erg->m_loc, erg->nat_loc, rotg->nat, erg->xc_center);
=======
            if ( (rotg->eType == erotgISOPF) || (rotg->eType == erotgPMPF) )
            {
                get_center_comm(cr, erg->x_loc_pbc, erg->m_loc, erg->nat_loc, rotg->nat, erg->xc_center);
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }

    } /* End of loop over rotation groups */

    /**************************************************************************/
    /* Done communicating, we can start to count cycles for the load balancing now ... */
    cycles_comp = gmx_cycles_read();


#ifdef TAKETIME
    t0 = MPI_Wtime();
#endif

<<<<<<< HEAD
    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg=rotg->enfrotgrp;
=======
    for (g = 0; g < rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg  = rotg->enfrotgrp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        bFlex = ISFLEX(rotg);
        bColl = ISCOLL(rotg);

        if (outstep_rot && MASTER(cr))
<<<<<<< HEAD
            fprintf(er->out_rot, "%12.4f", erg->degangle);
=======
        {
            fprintf(er->out_rot, "%12.4f", erg->degangle);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* Calculate angles and rotation matrices for potential fitting: */
        if ( (outstep_rot || outstep_slab) && (erotgFitPOT == rotg->eFittype) )
        {
            fit = erg->PotAngleFit;
            for (i = 0; i < rotg->PotAngle_nstep; i++)
            {
                calc_rotmat(rotg->vec, erg->degangle + fit->degangle[i], fit->rotmat[i]);

                /* Clear value from last step */
                erg->PotAngleFit->V[i] = 0.0;
            }
        }

        /* Clear values from last time step */
        erg->V        = 0.0;
        erg->torque_v = 0.0;
        erg->angle_v  = 0.0;
        erg->weight_v = 0.0;

<<<<<<< HEAD
        switch(rotg->eType)
=======
        switch (rotg->eType)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            case erotgISO:
            case erotgISOPF:
            case erotgPM:
            case erotgPMPF:
<<<<<<< HEAD
                do_fixed(rotg,x,box,t,step,outstep_rot,outstep_slab);
                break;
            case erotgRM:
                do_radial_motion(rotg,x,box,t,step,outstep_rot,outstep_slab);
                break;
            case erotgRMPF:
                do_radial_motion_pf(rotg,x,box,t,step,outstep_rot,outstep_slab);
                break;
            case erotgRM2:
            case erotgRM2PF:
                do_radial_motion2(rotg,x,box,t,step,outstep_rot,outstep_slab);
=======
                do_fixed(rotg, x, box, t, step, outstep_rot, outstep_slab);
                break;
            case erotgRM:
                do_radial_motion(rotg, x, box, t, step, outstep_rot, outstep_slab);
                break;
            case erotgRMPF:
                do_radial_motion_pf(rotg, x, box, t, step, outstep_rot, outstep_slab);
                break;
            case erotgRM2:
            case erotgRM2PF:
                do_radial_motion2(rotg, x, box, t, step, outstep_rot, outstep_slab);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                break;
            case erotgFLEXT:
            case erotgFLEX2T:
                /* Subtract the center of the rotation group from the collective positions array
                 * Also store the center in erg->xc_center since it needs to be subtracted
                 * in the low level routines from the local coordinates as well */
                get_center(erg->xc, erg->mc, rotg->nat, erg->xc_center);
                svmul(-1.0, erg->xc_center, transvec);
                translate_x(erg->xc, rotg->nat, transvec);
<<<<<<< HEAD
                do_flexible(MASTER(cr),er,rotg,g,x,box,t,step,outstep_rot,outstep_slab);
=======
                do_flexible(MASTER(cr), er, rotg, g, x, box, t, step, outstep_rot, outstep_slab);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                break;
            case erotgFLEX:
            case erotgFLEX2:
                /* Do NOT subtract the center of mass in the low level routines! */
                clear_rvec(erg->xc_center);
<<<<<<< HEAD
                do_flexible(MASTER(cr),er,rotg,g,x,box,t,step,outstep_rot,outstep_slab);
=======
                do_flexible(MASTER(cr), er, rotg, g, x, box, t, step, outstep_rot, outstep_slab);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                break;
            default:
                gmx_fatal(FARGS, "No such rotation potential.");
                break;
        }
    }

#ifdef TAKETIME
    if (MASTER(cr))
<<<<<<< HEAD
        fprintf(stderr, "%s calculation (step %d) took %g seconds.\n", RotStr, step, MPI_Wtime()-t0);
=======
    {
        fprintf(stderr, "%s calculation (step %d) took %g seconds.\n", RotStr, step, MPI_Wtime()-t0);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

    /* Stop the enforced rotation cycle counter and add the computation-only
     * cycles to the force cycles for load balancing */
    cycles_comp  = gmx_cycles_read() - cycles_comp;

    if (DOMAINDECOMP(cr) && wcycle)
<<<<<<< HEAD
        dd_cycles_add(cr->dd,cycles_comp,ddCyclF);
=======
    {
        dd_cycles_add(cr->dd, cycles_comp, ddCyclF);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}
