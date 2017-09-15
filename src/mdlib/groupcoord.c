<<<<<<< HEAD
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
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
 * Groningen Machine for Chemical Simulation
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "groupcoord.h"
#include "mpelogging.h"
#include "network.h"
#include "pbc.h"
#include "vec.h"
#include "smalloc.h"
#include "gmx_ga2la.h"

<<<<<<< HEAD
#define MIN(a,b) (((a)<(b))?(a):(b))



/* Select the indices of the group's atoms which are local and store them in 
=======
#define MIN(a, b) (((a) < (b)) ? (a) : (b))



/* Select the indices of the group's atoms which are local and store them in
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * anrs_loc[0..nr_loc]. The indices are saved in coll_ind[] for later reduction
 * in communicate_group_positions()
 */
extern void dd_make_local_group_indices(
<<<<<<< HEAD
        gmx_ga2la_t   ga2la,
        const int     nr,          /* IN:  Total number of atoms in the group */
        int           anrs[],      /* IN:  Global atom numbers of the groups atoms */
        int           *nr_loc,     /* OUT: Number of group atoms found locally */
        int           *anrs_loc[], /* OUT: Local atom numbers of the group  */
        int           *nalloc_loc, /* IN+OUT: Allocation size of anrs_loc */
        int           coll_ind[])  /* OUT (opt): Where is this position found in the collective array? */
{
    int  i,ii;
    int  localnr;       

    
    /* Loop over all the atom indices of the group to check
     * which ones are on the local node */
    localnr = 0;
    for(i=0; i<nr; i++)
    {
        if (ga2la_get_home(ga2la,anrs[i],&ii))
=======
        gmx_ga2la_t    ga2la,
        const int      nr,         /* IN:  Total number of atoms in the group */
        int            anrs[],     /* IN:  Global atom numbers of the groups atoms */
        int           *nr_loc,     /* OUT: Number of group atoms found locally */
        int           *anrs_loc[], /* OUT: Local atom numbers of the group  */
        int           *nalloc_loc, /* IN+OUT: Allocation size of anrs_loc */
        int            coll_ind[]) /* OUT (opt): Where is this position found in the collective array? */
{
    int  i, ii;
    int  localnr;


    /* Loop over all the atom indices of the group to check
     * which ones are on the local node */
    localnr = 0;
    for (i = 0; i < nr; i++)
    {
        if (ga2la_get_home(ga2la, anrs[i], &ii))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* The atom with this index is a home atom */
            if (localnr >= *nalloc_loc) /* Check whether memory suffices */
            {
                *nalloc_loc = over_alloc_dd(localnr+1);
                /* We never need more memory than the number of atoms in the group */
                *nalloc_loc = MIN(*nalloc_loc, nr);
<<<<<<< HEAD
                srenew(*anrs_loc,*nalloc_loc);
=======
                srenew(*anrs_loc, *nalloc_loc);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            /* Save the atoms index in the local atom numbers array */
            (*anrs_loc)[localnr] = ii;

            if (coll_ind != NULL)
            {
                /* Keep track of where this local atom belongs in the collective index array.
                 * This is needed when reducing the local arrays to a collective/global array
                 * in communicate_group_positions */
                coll_ind[localnr] = i;
            }

            /* add one to the local atom count */
            localnr++;
        }
    }
<<<<<<< HEAD
 
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Return the number of local atoms that were found */
    *nr_loc = localnr;
}


static void get_shifts_group(
<<<<<<< HEAD
        int    npbcdim, 
        matrix box,
        rvec   *xcoll,     /* IN:  Collective set of positions [0..nr] */
        int    nr,         /* IN:  Total number of atoms in the group */
        rvec   *xcoll_old, /* IN:  Positions from the last time step [0...nr] */
        ivec   *shifts)    /* OUT: Shifts for xcoll */
{
    int  i,m,d;
=======
        int     npbcdim,
        matrix  box,
        rvec   *xcoll,     /* IN:  Collective set of positions [0..nr] */
        int     nr,        /* IN:  Total number of atoms in the group */
        rvec   *xcoll_old, /* IN:  Positions from the last time step [0...nr] */
        ivec   *shifts)    /* OUT: Shifts for xcoll */
{
    int  i, m, d;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rvec dx;


    /* Get the shifts such that each atom is within closest
     * distance to its position at the last NS time step after shifting.
<<<<<<< HEAD
     * If we start with a whole group, and always keep track of 
     * shift changes, the group will stay whole this way */
    for (i=0; i < nr; i++)
        clear_ivec(shifts[i]);

    for (i=0; i<nr; i++)
    {
        /* The distance this atom moved since the last time step */
        /* If this is more than just a bit, it has changed its home pbc box */
        rvec_sub(xcoll[i],xcoll_old[i],dx);

        for(m=npbcdim-1; m>=0; m--)
        {
            while (dx[m] < -0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] += box[m][d];
=======
     * If we start with a whole group, and always keep track of
     * shift changes, the group will stay whole this way */
    for (i = 0; i < nr; i++)
    {
        clear_ivec(shifts[i]);
    }

    for (i = 0; i < nr; i++)
    {
        /* The distance this atom moved since the last time step */
        /* If this is more than just a bit, it has changed its home pbc box */
        rvec_sub(xcoll[i], xcoll_old[i], dx);

        for (m = npbcdim-1; m >= 0; m--)
        {
            while (dx[m] < -0.5*box[m][m])
            {
                for (d = 0; d < DIM; d++)
                {
                    dx[d] += box[m][d];
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                shifts[i][m]++;
            }
            while (dx[m] >= 0.5*box[m][m])
            {
<<<<<<< HEAD
                for(d=0; d<DIM; d++)
                    dx[d] -= box[m][d];
=======
                for (d = 0; d < DIM; d++)
                {
                    dx[d] -= box[m][d];
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                shifts[i][m]--;
            }
        }
    }
}


static void shift_positions_group(
<<<<<<< HEAD
        matrix box, 
        rvec   x[],      /* The positions [0..nr] */ 
        ivec   *is,      /* The shifts [0..nr] */ 
        int    nr)       /* The number of positions and shifts */
{
    int      i,tx,ty,tz;
=======
        matrix  box,
        rvec    x[],     /* The positions [0..nr] */
        ivec   *is,      /* The shifts [0..nr] */
        int     nr)      /* The number of positions and shifts */
{
    int      i, tx, ty, tz;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    GMX_MPE_LOG(ev_shift_start);

    /* Loop over the group's atoms */
<<<<<<< HEAD
    if(TRICLINIC(box)) 
    {
        for (i=0; i < nr; i++)
        {
            tx=is[i][XX];
            ty=is[i][YY];
            tz=is[i][ZZ];

            x[i][XX]=x[i][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
            x[i][YY]=x[i][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
            x[i][ZZ]=x[i][ZZ]+tz*box[ZZ][ZZ];
        }
    } else
    {
        for (i=0; i < nr; i++)
        {
            tx=is[i][XX];
            ty=is[i][YY];
            tz=is[i][ZZ];

            x[i][XX]=x[i][XX]+tx*box[XX][XX];
            x[i][YY]=x[i][YY]+ty*box[YY][YY];
            x[i][ZZ]=x[i][ZZ]+tz*box[ZZ][ZZ];
        }
    }    
=======
    if (TRICLINIC(box))
    {
        for (i = 0; i < nr; i++)
        {
            tx = is[i][XX];
            ty = is[i][YY];
            tz = is[i][ZZ];

            x[i][XX] = x[i][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
            x[i][YY] = x[i][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
            x[i][ZZ] = x[i][ZZ]+tz*box[ZZ][ZZ];
        }
    }
    else
    {
        for (i = 0; i < nr; i++)
        {
            tx = is[i][XX];
            ty = is[i][YY];
            tz = is[i][ZZ];

            x[i][XX] = x[i][XX]+tx*box[XX][XX];
            x[i][YY] = x[i][YY]+ty*box[YY][YY];
            x[i][ZZ] = x[i][ZZ]+tz*box[ZZ][ZZ];
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    GMX_MPE_LOG(ev_shift_finish);
}


<<<<<<< HEAD
/* Assemble the positions of the group such that every node has all of them. 
 * The atom indices are retrieved from anrs_loc[0..nr_loc] 
 * Note that coll_ind[i] = i is needed in the serial case */
extern void communicate_group_positions(
        t_commrec  *cr, 
        rvec       *xcoll,        /* OUT: Collective array of positions */
        ivec       *shifts,       /* IN+OUT: Collective array of shifts for xcoll */
        ivec       *extra_shifts, /* BUF: Extra shifts since last time step */
        const gmx_bool bNS,           /* IN:  NS step, the shifts have changed */
        rvec       *x_loc,        /* IN:  Local positions on this node */ 
        const int  nr,            /* IN:  Total number of atoms in the group */
        const int  nr_loc,        /* IN:  Local number of atoms in the group */
        int        *anrs_loc,     /* IN:  Local atom numbers */
        int        *coll_ind,     /* IN:  Collective index */
        rvec       *xcoll_old,    /* IN+OUT: Positions from the last time step, used to make group whole */
        matrix     box)
=======
/* Assemble the positions of the group such that every node has all of them.
 * The atom indices are retrieved from anrs_loc[0..nr_loc]
 * Note that coll_ind[i] = i is needed in the serial case */
extern void communicate_group_positions(
        t_commrec     *cr,
        rvec          *xcoll,        /* OUT: Collective array of positions */
        ivec          *shifts,       /* IN+OUT: Collective array of shifts for xcoll */
        ivec          *extra_shifts, /* BUF: Extra shifts since last time step */
        const gmx_bool bNS,          /* IN:  NS step, the shifts have changed */
        rvec          *x_loc,        /* IN:  Local positions on this node */
        const int      nr,           /* IN:  Total number of atoms in the group */
        const int      nr_loc,       /* IN:  Local number of atoms in the group */
        int           *anrs_loc,     /* IN:  Local atom numbers */
        int           *coll_ind,     /* IN:  Collective index */
        rvec          *xcoll_old,    /* IN+OUT: Positions from the last time step, used to make group whole */
        matrix         box)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    int i;


    GMX_MPE_LOG(ev_get_group_x_start);

    /* Zero out the groups' global position array */
    clear_rvecs(nr, xcoll);

<<<<<<< HEAD
    /* Put the local positions that this node has into the right place of 
     * the collective array. Note that in the serial case, coll_ind[i] = i */
    for (i=0; i<nr_loc; i++)
        copy_rvec(x_loc[anrs_loc[i]], xcoll[coll_ind[i]]);
=======
    /* Put the local positions that this node has into the right place of
     * the collective array. Note that in the serial case, coll_ind[i] = i */
    for (i = 0; i < nr_loc; i++)
    {
        copy_rvec(x_loc[anrs_loc[i]], xcoll[coll_ind[i]]);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (PAR(cr))
    {
        /* Add the arrays from all nodes together */
        gmx_sum(nr*3, xcoll[0], cr);
<<<<<<< HEAD

        /* To make the group whole, start with a whole group and each
         * step move the assembled positions at closest distance to the positions 
         * from the last step. First shift the positions with the saved shift 
         * vectors (these are 0 when this routine is called for the first time!) */
        shift_positions_group(box, xcoll, shifts, nr);
        
        /* Now check if some shifts changed since the last step.
         * This only needs to be done when the shifts are expected to have changed,
         * i.e. after neighboursearching */
        if (bNS) 
        {
            get_shifts_group(3, box, xcoll, nr, xcoll_old, extra_shifts);
            
            /* Shift with the additional shifts such that we get a whole group now */
            shift_positions_group(box, xcoll, extra_shifts, nr);
            
            /* Add the shift vectors together for the next time step */
            for (i=0; i<nr; i++)
            {
                shifts[i][XX] += extra_shifts[i][XX];
                shifts[i][YY] += extra_shifts[i][YY];
                shifts[i][ZZ] += extra_shifts[i][ZZ];
            }
            
            /* Store current correctly-shifted positions for comparison in the next NS time step */
            for (i=0; i<nr; i++)
                copy_rvec(xcoll[i],xcoll_old[i]);   
        }
    }
    
=======
    }
    /* To make the group whole, start with a whole group and each
     * step move the assembled positions at closest distance to the positions
     * from the last step. First shift the positions with the saved shift
     * vectors (these are 0 when this routine is called for the first time!) */
    shift_positions_group(box, xcoll, shifts, nr);

    /* Now check if some shifts changed since the last step.
     * This only needs to be done when the shifts are expected to have changed,
     * i.e. after neighboursearching */
    if (bNS)
    {
        get_shifts_group(3, box, xcoll, nr, xcoll_old, extra_shifts);

        /* Shift with the additional shifts such that we get a whole group now */
        shift_positions_group(box, xcoll, extra_shifts, nr);

        /* Add the shift vectors together for the next time step */
        for (i = 0; i < nr; i++)
        {
            shifts[i][XX] += extra_shifts[i][XX];
            shifts[i][YY] += extra_shifts[i][YY];
            shifts[i][ZZ] += extra_shifts[i][ZZ];
        }

        /* Store current correctly-shifted positions for comparison in the next NS time step */
        for (i = 0; i < nr; i++)
        {
            copy_rvec(xcoll[i], xcoll_old[i]);
        }
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    GMX_MPE_LOG(ev_get_group_x_finish);
}


/* Determine the (weighted) sum vector from positions x */
extern double get_sum_of_positions(rvec x[], real weight[], const int nat, dvec dsumvec)
{
<<<<<<< HEAD
    int i;
    rvec x_weighted;
=======
    int    i;
    rvec   x_weighted;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    double weight_sum = 0.0;


    /* Zero out the center */
    clear_dvec(dsumvec);

    /* Loop over all atoms and add their weighted position vectors */
    if (weight != NULL)
    {
<<<<<<< HEAD
        for (i=0; i<nat; i++)
=======
        for (i = 0; i < nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            weight_sum += weight[i];
            svmul(weight[i], x[i], x_weighted);
            dsumvec[XX] += x_weighted[XX];
            dsumvec[YY] += x_weighted[YY];
            dsumvec[ZZ] += x_weighted[ZZ];
        }
    }
    else
    {
<<<<<<< HEAD
        for (i=0; i<nat; i++)
=======
        for (i = 0; i < nat; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            dsumvec[XX] += x[i][XX];
            dsumvec[YY] += x[i][YY];
            dsumvec[ZZ] += x[i][ZZ];
        }
    }
    return weight_sum;
}


/* Determine center of structure from collective positions x */
extern void get_center(rvec x[], real weight[], const int nr, rvec rcenter)
{
    dvec   dcenter;
    double weight_sum, denom;

<<<<<<< HEAD
    
    weight_sum = get_sum_of_positions(x, weight, nr, dcenter);
    
    if (weight != NULL)
        denom = weight_sum; /* Divide by the sum of weight */
    else
        denom = nr;        /* Divide by the number of atoms */
        
    dsvmul(1.0/denom, dcenter, dcenter);
    
=======

    weight_sum = get_sum_of_positions(x, weight, nr, dcenter);

    if (weight != NULL)
    {
        denom = weight_sum; /* Divide by the sum of weight */
    }
    else
    {
        denom = nr;        /* Divide by the number of atoms */

    }
    dsvmul(1.0/denom, dcenter, dcenter);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rcenter[XX] = dcenter[XX];
    rcenter[YY] = dcenter[YY];
    rcenter[ZZ] = dcenter[ZZ];
}


/* Get the center from local positions that already have the correct
 * PBC representation */
extern void get_center_comm(
        t_commrec *cr,
<<<<<<< HEAD
        rvec x_loc[],       /* Local positions */
        real weight_loc[],  /* Local masses or other weights */
        int nr_loc,         /* Local number of atoms */
        int nr_group,       /* Total number of atoms of the group */ 
        rvec center)        /* Weighted center */
{
    double weight_sum, denom;
    dvec   dsumvec;
    double buf[4];    
    
    
    weight_sum = get_sum_of_positions(x_loc, weight_loc, nr_loc, dsumvec);
    
    /* Add the local contributions from all nodes. Put the sum vector and the 
=======
        rvec       x_loc[],      /* Local positions */
        real       weight_loc[], /* Local masses or other weights */
        int        nr_loc,       /* Local number of atoms */
        int        nr_group,     /* Total number of atoms of the group */
        rvec       center)       /* Weighted center */
{
    double weight_sum, denom;
    dvec   dsumvec;
    double buf[4];


    weight_sum = get_sum_of_positions(x_loc, weight_loc, nr_loc, dsumvec);

    /* Add the local contributions from all nodes. Put the sum vector and the
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
     * weight in a buffer array so that we get along with a single communication
     * call. */
    if (PAR(cr))
    {
        buf[0] = dsumvec[XX];
        buf[1] = dsumvec[YY];
        buf[2] = dsumvec[ZZ];
        buf[3] = weight_sum;
<<<<<<< HEAD
        
        /* Communicate buffer */
        gmx_sumd(4, buf, cr);
        
=======

        /* Communicate buffer */
        gmx_sumd(4, buf, cr);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        dsumvec[XX] = buf[0];
        dsumvec[YY] = buf[1];
        dsumvec[ZZ] = buf[2];
        weight_sum  = buf[3];
    }
<<<<<<< HEAD
    
    if (weight_loc != NULL)
        denom = 1.0/weight_sum; /* Divide by the sum of weight to get center of mass e.g. */
    else
        denom = 1.0/nr_group;   /* Divide by the number of atoms to get the geometrical center */
        
=======

    if (weight_loc != NULL)
    {
        denom = 1.0/weight_sum; /* Divide by the sum of weight to get center of mass e.g. */
    }
    else
    {
        denom = 1.0/nr_group;   /* Divide by the number of atoms to get the geometrical center */

    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    center[XX] = dsumvec[XX]*denom;
    center[YY] = dsumvec[YY]*denom;
    center[ZZ] = dsumvec[ZZ]*denom;
}


/* Translate x with transvec */
extern void translate_x(rvec x[], const int nr, const rvec transvec)
{
    int i;
<<<<<<< HEAD
    
    
    for (i=0; i<nr; i++)
        rvec_inc(x[i], transvec);
=======


    for (i = 0; i < nr; i++)
    {
        rvec_inc(x[i], transvec);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


extern void rotate_x(rvec x[], const int nr, matrix rmat)
{
<<<<<<< HEAD
    int i,j,k;
    rvec x_old;

    
    /* Apply the rotation matrix */
    for (i=0; i<nr; i++)
    {
        for (j=0; j<3; j++)
            x_old[j] = x[i][j];
        for (j=0; j<3; j++)
        {
            x[i][j] = 0;
            for (k=0; k<3; k++)
                x[i][j] += rmat[k][j]*x_old[k];
        }
    }
}

=======
    int  i, j, k;
    rvec x_old;


    /* Apply the rotation matrix */
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x_old[j] = x[i][j];
        }
        for (j = 0; j < 3; j++)
        {
            x[i][j] = 0;
            for (k = 0; k < 3; k++)
            {
                x[i][j] += rmat[k][j]*x_old[k];
            }
        }
    }
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
