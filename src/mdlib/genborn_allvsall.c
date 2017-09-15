/*
<<<<<<< HEAD
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
=======
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
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

#include <math.h>
#include "types/simple.h"

#include "vec.h"
#include "smalloc.h"

#include "partdec.h"
#include "network.h"
#include "physics.h"
#include "genborn.h"
#include "genborn_allvsall.h"


<<<<<<< HEAD
typedef struct 
{
    int *      jindex_gb;
    int **     exclusion_mask_gb;
} 
gmx_allvsallgb2_data_t;

static int 
calc_maxoffset(int i,int natoms)
{
    int maxoffset;
    
=======
typedef struct
{
    int *      jindex_gb;
    int **     exclusion_mask_gb;
}
gmx_allvsallgb2_data_t;

static int
calc_maxoffset(int i, int natoms)
{
    int maxoffset;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if ((natoms % 2) == 1)
    {
        /* Odd number of atoms, easy */
        maxoffset = natoms/2;
    }
    else if ((natoms % 4) == 0)
    {
        /* Multiple of four is hard */
        if (i < natoms/2)
        {
            if ((i % 2) == 0)
            {
<<<<<<< HEAD
                maxoffset=natoms/2;
            }
            else
            {
                maxoffset=natoms/2-1;
=======
                maxoffset = natoms/2;
            }
            else
            {
                maxoffset = natoms/2-1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
        else
        {
            if ((i % 2) == 1)
            {
<<<<<<< HEAD
                maxoffset=natoms/2;
            }
            else
            {
                maxoffset=natoms/2-1;
=======
                maxoffset = natoms/2;
            }
            else
            {
                maxoffset = natoms/2-1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
    else
    {
        /* natoms/2 = odd */
        if ((i % 2) == 0)
        {
<<<<<<< HEAD
            maxoffset=natoms/2;
        }
        else
        {
            maxoffset=natoms/2-1;
        }
    }
    
=======
            maxoffset = natoms/2;
        }
        else
        {
            maxoffset = natoms/2-1;
        }
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return maxoffset;
}

static void
<<<<<<< HEAD
setup_gb_exclusions_and_indices(gmx_allvsallgb2_data_t *   aadata,
                                t_ilist *                  ilist,
                                int                        natoms,
=======
setup_gb_exclusions_and_indices(gmx_allvsallgb2_data_t     *   aadata,
                                t_ilist     *                  ilist,
                                int                            natoms,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                gmx_bool                       bInclude12,
                                gmx_bool                       bInclude13,
                                gmx_bool                       bInclude14)
{
<<<<<<< HEAD
    int i,j,k,tp;
    int a1,a2;
    int nj0,nj1;
    int max_offset;
    int max_excl_offset;
    int nj;
    
    /* This routine can appear to be a bit complex, but it is mostly book-keeping.
     * To enable the fast all-vs-all kernel we need to be able to stream through all coordinates
     * whether they should interact or not. 
=======
    int i, j, k, tp;
    int a1, a2;
    int nj0, nj1;
    int max_offset;
    int max_excl_offset;
    int nj;

    /* This routine can appear to be a bit complex, but it is mostly book-keeping.
     * To enable the fast all-vs-all kernel we need to be able to stream through all coordinates
     * whether they should interact or not.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
     *
     * To avoid looping over the exclusions, we create a simple mask that is 1 if the interaction
     * should be present, otherwise 0. Since exclusions typically only occur when i & j are close,
     * we create a jindex array with three elements per i atom: the starting point, the point to
     * which we need to check exclusions, and the end point.
     * This way we only have to allocate a short exclusion mask per i atom.
     */
<<<<<<< HEAD
    
    /* Allocate memory for jindex arrays */
    snew(aadata->jindex_gb,3*natoms);
    
    /* Pointer to lists with exclusion masks */
    snew(aadata->exclusion_mask_gb,natoms);
    
    for(i=0;i<natoms;i++)
    {
        /* Start */
        aadata->jindex_gb[3*i]       = i+1;        
        max_offset = calc_maxoffset(i,natoms);

        /* first check the max range of atoms to EXCLUDE */
        max_excl_offset = 0;
        if(!bInclude12)
        {
            for(j=0;j<ilist[F_GB12].nr;j+=3)
            {
                a1 = ilist[F_GB12].iatoms[j+1];
                a2 = ilist[F_GB12].iatoms[j+2];
                
                if(a1==i)
                {
                    k = a2-a1;
                }
                else if(a2==i)
                {
                    k = a1+natoms-a2;
                }
                else 
                {
                    continue;
                }
                if(k>0 && k<=max_offset)
=======

    /* Allocate memory for jindex arrays */
    snew(aadata->jindex_gb, 3*natoms);

    /* Pointer to lists with exclusion masks */
    snew(aadata->exclusion_mask_gb, natoms);

    for (i = 0; i < natoms; i++)
    {
        /* Start */
        aadata->jindex_gb[3*i]       = i+1;
        max_offset                   = calc_maxoffset(i, natoms);

        /* first check the max range of atoms to EXCLUDE */
        max_excl_offset = 0;
        if (!bInclude12)
        {
            for (j = 0; j < ilist[F_GB12].nr; j += 3)
            {
                a1 = ilist[F_GB12].iatoms[j+1];
                a2 = ilist[F_GB12].iatoms[j+2];

                if (a1 == i)
                {
                    k = a2-a1;
                }
                else if (a2 == i)
                {
                    k = a1+natoms-a2;
                }
                else
                {
                    continue;
                }
                if (k > 0 && k <= max_offset)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
                }
            }
        }
<<<<<<< HEAD
        if(!bInclude13)
        {
            for(j=0;j<ilist[F_GB13].nr;j+=3)
            {
                a1 = ilist[F_GB13].iatoms[j+1];
                a2 = ilist[F_GB13].iatoms[j+2];
                
                
                if(a1==i)
                {
                    k = a2-a1;
                }
                else if(a2==i)
                {
                    k = a1+natoms-a2;
                }
                else 
                {
                    continue;
                }
                if(k>0 && k<=max_offset)
=======
        if (!bInclude13)
        {
            for (j = 0; j < ilist[F_GB13].nr; j += 3)
            {
                a1 = ilist[F_GB13].iatoms[j+1];
                a2 = ilist[F_GB13].iatoms[j+2];


                if (a1 == i)
                {
                    k = a2-a1;
                }
                else if (a2 == i)
                {
                    k = a1+natoms-a2;
                }
                else
                {
                    continue;
                }
                if (k > 0 && k <= max_offset)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
                }
            }
        }
<<<<<<< HEAD
        if(!bInclude14)
        {
            for(j=0;j<ilist[F_GB14].nr;j+=3)
            {
                a1 = ilist[F_GB14].iatoms[j+1];
                a2 = ilist[F_GB14].iatoms[j+2];
                
                
                if(a1==i)
                {
                    k = a2-a1;
                }
                else if(a2==i)
                {
                    k = a1+natoms-a2;
                }
                else 
                {
                    continue;
                }
                if(k>0 && k<=max_offset)
=======
        if (!bInclude14)
        {
            for (j = 0; j < ilist[F_GB14].nr; j += 3)
            {
                a1 = ilist[F_GB14].iatoms[j+1];
                a2 = ilist[F_GB14].iatoms[j+2];


                if (a1 == i)
                {
                    k = a2-a1;
                }
                else if (a2 == i)
                {
                    k = a1+natoms-a2;
                }
                else
                {
                    continue;
                }
                if (k > 0 && k <= max_offset)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
                }
            }
        }
        max_excl_offset = (max_offset < max_excl_offset) ? max_offset : max_excl_offset;

<<<<<<< HEAD
        aadata->jindex_gb[3*i+1] = i+1+max_excl_offset;        
        
        snew(aadata->exclusion_mask_gb[i],max_excl_offset);
        
        /* Include everything by default */
        for(j=0;j<max_excl_offset;j++)
=======
        aadata->jindex_gb[3*i+1] = i+1+max_excl_offset;

        snew(aadata->exclusion_mask_gb[i], max_excl_offset);

        /* Include everything by default */
        for (j = 0; j < max_excl_offset; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* Use all-ones to mark interactions that should be present, compatible with SSE */
            aadata->exclusion_mask_gb[i][j] = 0xFFFFFFFF;
        }
        /* Go through exclusions again */
<<<<<<< HEAD
        if(!bInclude12)
        {
            for(j=0;j<ilist[F_GB12].nr;j+=3)
            {
                a1 = ilist[F_GB12].iatoms[j+1];
                a2 = ilist[F_GB12].iatoms[j+2];
                
                if(a1==i)
                {
                    k = a2-a1;
                }
                else if(a2==i)
                {
                    k = a1+natoms-a2;
                }
                else 
                {
                    continue;
                }
                if(k>0 && k<=max_offset)
=======
        if (!bInclude12)
        {
            for (j = 0; j < ilist[F_GB12].nr; j += 3)
            {
                a1 = ilist[F_GB12].iatoms[j+1];
                a2 = ilist[F_GB12].iatoms[j+2];

                if (a1 == i)
                {
                    k = a2-a1;
                }
                else if (a2 == i)
                {
                    k = a1+natoms-a2;
                }
                else
                {
                    continue;
                }
                if (k > 0 && k <= max_offset)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    aadata->exclusion_mask_gb[i][k-1] = 0;
                }
            }
        }
<<<<<<< HEAD
        if(!bInclude13)
        {
            for(j=0;j<ilist[F_GB13].nr;j+=3)
            {
                a1 = ilist[F_GB13].iatoms[j+1];
                a2 = ilist[F_GB13].iatoms[j+2];
               
                if(a1==i)
                {
                    k = a2-a1;
                }
                else if(a2==i)
                {
                    k = a1+natoms-a2;
                }
                else 
                {
                    continue;
                }
                if(k>0 && k<=max_offset)
=======
        if (!bInclude13)
        {
            for (j = 0; j < ilist[F_GB13].nr; j += 3)
            {
                a1 = ilist[F_GB13].iatoms[j+1];
                a2 = ilist[F_GB13].iatoms[j+2];

                if (a1 == i)
                {
                    k = a2-a1;
                }
                else if (a2 == i)
                {
                    k = a1+natoms-a2;
                }
                else
                {
                    continue;
                }
                if (k > 0 && k <= max_offset)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    aadata->exclusion_mask_gb[i][k-1] = 0;
                }
            }
        }
<<<<<<< HEAD
        if(!bInclude14)
        {
            for(j=0;j<ilist[F_GB14].nr;j+=3)
=======
        if (!bInclude14)
        {
            for (j = 0; j < ilist[F_GB14].nr; j += 3)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                a1 = ilist[F_GB14].iatoms[j+1];
                a2 = ilist[F_GB14].iatoms[j+2];

<<<<<<< HEAD
                if(a1==i)
                {
                    k = a2-a1;
                }
                else if(a2==i)
                {
                    k = a1+natoms-a2;
                }
                else 
                {
                    continue;
                }
                if(k>0 && k<=max_offset)
=======
                if (a1 == i)
                {
                    k = a2-a1;
                }
                else if (a2 == i)
                {
                    k = a1+natoms-a2;
                }
                else
                {
                    continue;
                }
                if (k > 0 && k <= max_offset)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    aadata->exclusion_mask_gb[i][k-1] = 0;
                }
            }
        }
<<<<<<< HEAD
        
        /* End */
        
        /* End */
        aadata->jindex_gb[3*i+2] = i+1+max_offset;        
=======

        /* End */

        /* End */
        aadata->jindex_gb[3*i+2] = i+1+max_offset;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}


static void
<<<<<<< HEAD
genborn_allvsall_setup(gmx_allvsallgb2_data_t **  p_aadata,
                       t_ilist *                  ilist,
                       int                        natoms,
=======
genborn_allvsall_setup(gmx_allvsallgb2_data_t     **  p_aadata,
                       t_ilist     *                  ilist,
                       int                            natoms,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                       gmx_bool                       bInclude12,
                       gmx_bool                       bInclude13,
                       gmx_bool                       bInclude14)
{
<<<<<<< HEAD
	int i,j,idx;
	gmx_allvsallgb2_data_t *aadata;
    real *p;
    
	snew(aadata,1);
	*p_aadata = aadata;
        
    setup_gb_exclusions_and_indices(aadata,ilist,natoms,bInclude12,bInclude13,bInclude14);
=======
    int                     i, j, idx;
    gmx_allvsallgb2_data_t *aadata;
    real                   *p;

    snew(aadata, 1);
    *p_aadata = aadata;

    setup_gb_exclusions_and_indices(aadata, ilist, natoms, bInclude12, bInclude13, bInclude14);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}



int
genborn_allvsall_calc_still_radii(t_forcerec *           fr,
                                  t_mdatoms *            mdatoms,
                                  gmx_genborn_t *        born,
                                  gmx_localtop_t *       top,
                                  real *                 x,
                                  t_commrec *            cr,
                                  void *                 work)
{
<<<<<<< HEAD
	gmx_allvsallgb2_data_t *aadata;
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2;
	int        i,j,k,n;
    int *      mask;
    
    real       ix,iy,iz;
    real       jx,jy,jz;
    real       dx,dy,dz;
    real       rsq,rinv;
    real       gpi,rai,vai;
    real       prod_ai;
    real       irsq,idr4,idr6;
    real       raj,rvdw,ratio;
    real       vaj,ccf,dccf,theta,cosq;
    real       term,prod,icf4,icf6,gpi2,factor,sinq;
    
    natoms              = mdatoms->nr;
	ni0                 = mdatoms->start;
	ni1                 = mdatoms->start+mdatoms->homenr;
    factor  = 0.5*ONE_4PI_EPS0;
    n = 0;
    
    aadata = *((gmx_allvsallgb2_data_t **)work);
    
	if(aadata==NULL)
	{
		genborn_allvsall_setup(&aadata,top->idef.il,mdatoms->nr,
                               FALSE,FALSE,TRUE);
        *((gmx_allvsallgb2_data_t **)work) = aadata;
	}
    
    
    for(i=0;i<born->nr;i++)
    {
        born->gpol_still_work[i] = 0;
    }
    
    
	for(i=ni0; i<ni1; i++)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */

		/* Load i atom data */
        ix                = x[3*i];
        iy                = x[3*i+1];
        iz                = x[3*i+2];
        
        gpi               = 0.0;

        rai     = top->atomtypes.gb_radius[mdatoms->typeA[i]];
		vai     = born->vsolv[i];
		prod_ai = STILL_P4*vai;
        
		/* Load limits for loop over neighbors */
		nj0              = aadata->jindex_gb[3*i];
		nj1              = aadata->jindex_gb[3*i+1];
		nj2              = aadata->jindex_gb[3*i+2];
        
        mask             = aadata->exclusion_mask_gb[i];

        /* Prologue part, including exclusion mask */
        for(j=nj0; j<nj1; j++,mask++)
        {          
            if(*mask!=0)
=======
    gmx_allvsallgb2_data_t *aadata;
    int                     natoms;
    int                     ni0, ni1;
    int                     nj0, nj1, nj2;
    int                     i, j, k, n;
    int              *      mask;

    real                    ix, iy, iz;
    real                    jx, jy, jz;
    real                    dx, dy, dz;
    real                    rsq, rinv;
    real                    gpi, rai, vai;
    real                    prod_ai;
    real                    irsq, idr4, idr6;
    real                    raj, rvdw, ratio;
    real                    vaj, ccf, dccf, theta, cosq;
    real                    term, prod, icf4, icf6, gpi2, factor, sinq;

    natoms              = mdatoms->nr;
    ni0                 = mdatoms->start;
    ni1                 = mdatoms->start+mdatoms->homenr;
    factor              = 0.5*ONE_4PI_EPS0;
    n                   = 0;

    aadata = *((gmx_allvsallgb2_data_t **)work);

    if (aadata == NULL)
    {
        genborn_allvsall_setup(&aadata, top->idef.il, mdatoms->nr,
                               FALSE, FALSE, TRUE);
        *((gmx_allvsallgb2_data_t **)work) = aadata;
    }


    for (i = 0; i < born->nr; i++)
    {
        born->gpol_still_work[i] = 0;
    }


    for (i = ni0; i < ni1; i++)
    {
        /* We assume shifts are NOT used for all-vs-all interactions */

        /* Load i atom data */
        ix                = x[3*i];
        iy                = x[3*i+1];
        iz                = x[3*i+2];

        gpi               = 0.0;

        rai     = top->atomtypes.gb_radius[mdatoms->typeA[i]];
        vai     = born->vsolv[i];
        prod_ai = STILL_P4*vai;

        /* Load limits for loop over neighbors */
        nj0              = aadata->jindex_gb[3*i];
        nj1              = aadata->jindex_gb[3*i+1];
        nj2              = aadata->jindex_gb[3*i+2];

        mask             = aadata->exclusion_mask_gb[i];

        /* Prologue part, including exclusion mask */
        for (j = nj0; j < nj1; j++, mask++)
        {
            if (*mask != 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                k = j%natoms;

                /* load j atom coordinates */
                jx                = x[3*k];
                jy                = x[3*k+1];
                jz                = x[3*k+2];
<<<<<<< HEAD
                
                /* Calculate distance */
                dx                = ix - jx;      
                dy                = iy - jy;      
                dz                = iz - jz;      
                rsq               = dx*dx+dy*dy+dz*dz;
                
                /* Calculate 1/r and 1/r2 */
                rinv              = gmx_invsqrt(rsq);
                irsq  = rinv*rinv;
                idr4  = irsq*irsq;
                idr6  = idr4*irsq;

                raj = top->atomtypes.gb_radius[mdatoms->typeA[k]];
                
                rvdw  = rai + raj;
                
                ratio = rsq / (rvdw * rvdw);
                vaj   = born->vsolv[k];
                
                
                if(ratio>STILL_P5INV) 
                {
                    ccf=1.0;
                    dccf=0.0;
=======

                /* Calculate distance */
                dx                = ix - jx;
                dy                = iy - jy;
                dz                = iz - jz;
                rsq               = dx*dx+dy*dy+dz*dz;

                /* Calculate 1/r and 1/r2 */
                rinv              = gmx_invsqrt(rsq);
                irsq              = rinv*rinv;
                idr4              = irsq*irsq;
                idr6              = idr4*irsq;

                raj = top->atomtypes.gb_radius[mdatoms->typeA[k]];

                rvdw  = rai + raj;

                ratio = rsq / (rvdw * rvdw);
                vaj   = born->vsolv[k];


                if (ratio > STILL_P5INV)
                {
                    ccf  = 1.0;
                    dccf = 0.0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
                else
                {
                    theta = ratio*STILL_PIP5;
                    cosq  = cos(theta);
                    term  = 0.5*(1.0-cosq);
                    ccf   = term*term;
                    sinq  = 1.0 - cosq*cosq;
                    dccf  = 2.0*term*sinq*gmx_invsqrt(sinq)*theta;
                }
<<<<<<< HEAD
                
                prod          = STILL_P4*vaj;
                icf4          = ccf*idr4;
                icf6          = (4*ccf-dccf)*idr6;
                
                born->gpol_still_work[k] += prod_ai*icf4;
                gpi             = gpi+prod*icf4;
                
                /* Save ai->aj and aj->ai chain rule terms */
                fr->dadx[n++]   = prod*icf6;
                fr->dadx[n++]   = prod_ai*icf6;   
                
=======

                prod          = STILL_P4*vaj;
                icf4          = ccf*idr4;
                icf6          = (4*ccf-dccf)*idr6;

                born->gpol_still_work[k] += prod_ai*icf4;
                gpi                       = gpi+prod*icf4;

                /* Save ai->aj and aj->ai chain rule terms */
                fr->dadx[n++]   = prod*icf6;
                fr->dadx[n++]   = prod_ai*icf6;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                /* 27 flops, plus one cos(x) - estimate at 20 flops  => 47 */

            }
        }
<<<<<<< HEAD
        
        /* Main part, no exclusions */
        for(j=nj1; j<nj2; j++)
        {       
=======

        /* Main part, no exclusions */
        for (j = nj1; j < nj2; j++)
        {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            k = j%natoms;

            /* load j atom coordinates */
            jx                = x[3*k];
            jy                = x[3*k+1];
            jz                = x[3*k+2];
<<<<<<< HEAD
            
            /* Calculate distance */
            dx                = ix - jx;      
            dy                = iy - jy;      
            dz                = iz - jz;      
            rsq               = dx*dx+dy*dy+dz*dz;
            
            /* Calculate 1/r and 1/r2 */
            rinv              = gmx_invsqrt(rsq);
            irsq  = rinv*rinv;
            idr4  = irsq*irsq;
            idr6  = idr4*irsq;
            
            raj = top->atomtypes.gb_radius[mdatoms->typeA[k]];
            
            rvdw  = rai + raj;
            
            ratio = rsq / (rvdw * rvdw);
            vaj   = born->vsolv[k];
            
            if(ratio>STILL_P5INV) 
            {
                ccf=1.0;
                dccf=0.0;
=======

            /* Calculate distance */
            dx                = ix - jx;
            dy                = iy - jy;
            dz                = iz - jz;
            rsq               = dx*dx+dy*dy+dz*dz;

            /* Calculate 1/r and 1/r2 */
            rinv              = gmx_invsqrt(rsq);
            irsq              = rinv*rinv;
            idr4              = irsq*irsq;
            idr6              = idr4*irsq;

            raj = top->atomtypes.gb_radius[mdatoms->typeA[k]];

            rvdw  = rai + raj;

            ratio = rsq / (rvdw * rvdw);
            vaj   = born->vsolv[k];

            if (ratio > STILL_P5INV)
            {
                ccf  = 1.0;
                dccf = 0.0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            else
            {
                theta = ratio*STILL_PIP5;
                cosq  = cos(theta);
                term  = 0.5*(1.0-cosq);
                ccf   = term*term;
                sinq  = 1.0 - cosq*cosq;
                dccf  = 2.0*term*sinq*gmx_invsqrt(sinq)*theta;
            }
<<<<<<< HEAD
            
            prod          = STILL_P4*vaj;
            icf4          = ccf*idr4;
            icf6          = (4*ccf-dccf)*idr6;
            
            born->gpol_still_work[k] += prod_ai*icf4;
            gpi             = gpi+prod*icf4;
            
=======

            prod          = STILL_P4*vaj;
            icf4          = ccf*idr4;
            icf6          = (4*ccf-dccf)*idr6;

            born->gpol_still_work[k] += prod_ai*icf4;
            gpi                       = gpi+prod*icf4;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Save ai->aj and aj->ai chain rule terms */
            fr->dadx[n++]   = prod*icf6;
            fr->dadx[n++]   = prod_ai*icf6;
        }
<<<<<<< HEAD
        born->gpol_still_work[i]+=gpi;
	}    
    
    /* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_still_work, cr);
	}
	
	/* Calculate the radii */
	for(i=0;i<natoms;i++)
	{
		if(born->use[i] != 0)
		{
			gpi  = born->gpol[i]+born->gpol_still_work[i];
			gpi2 = gpi * gpi;
			born->bRad[i]   = factor*gmx_invsqrt(gpi2);
			fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
		}
	}
	
	return 0;
=======
        born->gpol_still_work[i] += gpi;
    }

    /* Parallel summations */
    if (PARTDECOMP(cr))
    {
        gmx_sum(natoms, born->gpol_still_work, cr);
    }

    /* Calculate the radii */
    for (i = 0; i < natoms; i++)
    {
        if (born->use[i] != 0)
        {
            gpi             = born->gpol[i]+born->gpol_still_work[i];
            gpi2            = gpi * gpi;
            born->bRad[i]   = factor*gmx_invsqrt(gpi2);
            fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
        }
    }

    return 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}



int
genborn_allvsall_calc_hct_obc_radii(t_forcerec *           fr,
                                    t_mdatoms *            mdatoms,
                                    gmx_genborn_t *        born,
                                    int                    gb_algorithm,
                                    gmx_localtop_t *       top,
                                    real *                 x,
                                    t_commrec *            cr,
                                    void *                 work)
{
<<<<<<< HEAD
	gmx_allvsallgb2_data_t *aadata;
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2;
	int        i,j,k,n;
    int *      mask;
    
    real       ix,iy,iz;
    real       jx,jy,jz;
    real       dx,dy,dz;
    real       rsq,rinv;
    real       prod,raj;
    real       rai,doffset,rai_inv,rai_inv2,sk_ai,sk2_ai,sum_ai;
    real       dr,sk,lij,dlij,lij2,lij3,uij2,uij3,diff2,uij,log_term;
    real       lij_inv,sk2,sk2_rinv,tmp,t1,t2,t3,raj_inv,sum_ai2,sum_ai3,tsum;
    real       tchain;
    real       dadxi,dadxj;
    real       rad,min_rad;
    
    natoms              = mdatoms->nr;
	ni0                 = mdatoms->start;
	ni1                 = mdatoms->start+mdatoms->homenr;

    n = 0;
    prod = 0;
    raj = 0;
    doffset = born->gb_doffset;

    aadata = *((gmx_allvsallgb2_data_t **)work);
    
	if(aadata==NULL)
	{
		genborn_allvsall_setup(&aadata,top->idef.il,mdatoms->nr,
                               TRUE,TRUE,TRUE);
        *((gmx_allvsallgb2_data_t **)work) = aadata;
	}
    
    for(i=0;i<born->nr;i++)
    {
        born->gpol_hct_work[i] = 0;
    }
    
	for(i=ni0; i<ni1; i++)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */
		
		/* Load i atom data */
        ix                = x[3*i];
        iy                = x[3*i+1];
        iz                = x[3*i+2];
        
		rai      = top->atomtypes.gb_radius[mdatoms->typeA[i]]-doffset;
		rai_inv  = 1.0/rai;
		
		sk_ai    = born->param[i];
		sk2_ai   = sk_ai*sk_ai;

        sum_ai   = 0;
        
		/* Load limits for loop over neighbors */
		nj0              = aadata->jindex_gb[3*i];
		nj1              = aadata->jindex_gb[3*i+1];
		nj2              = aadata->jindex_gb[3*i+2];
        
        mask             = aadata->exclusion_mask_gb[i];
               
        /* Prologue part, including exclusion mask */
        for(j=nj0; j<nj1; j++,mask++)
        {          
            if(*mask!=0)
            {
                k = j%natoms;
                
=======
    gmx_allvsallgb2_data_t *aadata;
    int                     natoms;
    int                     ni0, ni1;
    int                     nj0, nj1, nj2;
    int                     i, j, k, n;
    int              *      mask;

    real                    ix, iy, iz;
    real                    jx, jy, jz;
    real                    dx, dy, dz;
    real                    rsq, rinv;
    real                    prod, raj;
    real                    rai, doffset, rai_inv, rai_inv2, sk_ai, sk2_ai, sum_ai;
    real                    dr, sk, lij, dlij, lij2, lij3, uij2, uij3, diff2, uij, log_term;
    real                    lij_inv, sk2, sk2_rinv, tmp, t1, t2, t3, raj_inv, sum_ai2, sum_ai3, tsum;
    real                    tchain;
    real                    dadxi, dadxj;
    real                    rad, min_rad;

    natoms              = mdatoms->nr;
    ni0                 = mdatoms->start;
    ni1                 = mdatoms->start+mdatoms->homenr;

    n       = 0;
    prod    = 0;
    raj     = 0;
    doffset = born->gb_doffset;

    aadata = *((gmx_allvsallgb2_data_t **)work);

    if (aadata == NULL)
    {
        genborn_allvsall_setup(&aadata, top->idef.il, mdatoms->nr,
                               TRUE, TRUE, TRUE);
        *((gmx_allvsallgb2_data_t **)work) = aadata;
    }

    for (i = 0; i < born->nr; i++)
    {
        born->gpol_hct_work[i] = 0;
    }

    for (i = ni0; i < ni1; i++)
    {
        /* We assume shifts are NOT used for all-vs-all interactions */

        /* Load i atom data */
        ix                = x[3*i];
        iy                = x[3*i+1];
        iz                = x[3*i+2];

        rai      = top->atomtypes.gb_radius[mdatoms->typeA[i]]-doffset;
        rai_inv  = 1.0/rai;

        sk_ai    = born->param[i];
        sk2_ai   = sk_ai*sk_ai;

        sum_ai   = 0;

        /* Load limits for loop over neighbors */
        nj0              = aadata->jindex_gb[3*i];
        nj1              = aadata->jindex_gb[3*i+1];
        nj2              = aadata->jindex_gb[3*i+2];

        mask             = aadata->exclusion_mask_gb[i];

        /* Prologue part, including exclusion mask */
        for (j = nj0; j < nj1; j++, mask++)
        {
            if (*mask != 0)
            {
                k = j%natoms;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                /* load j atom coordinates */
                jx                = x[3*k];
                jy                = x[3*k+1];
                jz                = x[3*k+2];
<<<<<<< HEAD
                
                /* Calculate distance */
                dx                = ix - jx;      
                dy                = iy - jy;      
                dz                = iz - jz;      
                rsq               = dx*dx+dy*dy+dz*dz;
                
                /* Calculate 1/r and 1/r2 */
                rinv              = gmx_invsqrt(rsq);
                dr                = rsq*rinv;
                
                /* sk is precalculated in init_gb() */
                sk    = born->param[k];
                raj   = top->atomtypes.gb_radius[mdatoms->typeA[k]]-doffset; 
                
                /* aj -> ai interaction */
                               
                
                if(rai < dr+sk)
                {
                    lij       = 1.0/(dr-sk); 
                    dlij      = 1.0; 
                    
                    if(rai>dr-sk)
                    {
                        lij  = rai_inv; 
                        dlij = 0.0;
                    }
                                         
                    uij      = 1.0/(dr+sk);  
=======

                /* Calculate distance */
                dx                = ix - jx;
                dy                = iy - jy;
                dz                = iz - jz;
                rsq               = dx*dx+dy*dy+dz*dz;

                /* Calculate 1/r and 1/r2 */
                rinv              = gmx_invsqrt(rsq);
                dr                = rsq*rinv;

                /* sk is precalculated in init_gb() */
                sk    = born->param[k];
                raj   = top->atomtypes.gb_radius[mdatoms->typeA[k]]-doffset;

                /* aj -> ai interaction */


                if (rai < dr+sk)
                {
                    lij       = 1.0/(dr-sk);
                    dlij      = 1.0;

                    if (rai > dr-sk)
                    {
                        lij  = rai_inv;
                        dlij = 0.0;
                    }

                    uij      = 1.0/(dr+sk);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    lij2     = lij  * lij;
                    lij3     = lij2 * lij;
                    uij2     = uij  * uij;
                    uij3     = uij2 * uij;

<<<<<<< HEAD
                    diff2    = uij2-lij2; 
                    
                    lij_inv  = gmx_invsqrt(lij2); 
                    sk2      = sk*sk;
                    sk2_rinv = sk2*rinv;	
                    prod     = 0.25*sk2_rinv;
                    
                    log_term = log(uij*lij_inv); 
                    /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                    tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
                    
                    if(rai < sk-dr)
                    {
                        tmp = tmp + 2.0 * (rai_inv-lij); 
                    }
                    
                    t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);   
                    t2      = -0.5*uij2 - prod*uij3 + 0.25*(uij*rinv+uij3*dr); 
                    
                    t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
					
                    dadxi = (dlij*t1+t2+t3)*rinv; 
                    
=======
                    diff2    = uij2-lij2;

                    lij_inv  = gmx_invsqrt(lij2);
                    sk2      = sk*sk;
                    sk2_rinv = sk2*rinv;
                    prod     = 0.25*sk2_rinv;

                    log_term = log(uij*lij_inv);
                    /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                    tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);

                    if (rai < sk-dr)
                    {
                        tmp = tmp + 2.0 * (rai_inv-lij);
                    }

                    t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                    t2      = -0.5*uij2 - prod*uij3 + 0.25*(uij*rinv+uij3*dr);

                    t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                    dadxi = (dlij*t1+t2+t3)*rinv;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    sum_ai += 0.5*tmp;
                }
                else
                {
                    dadxi = 0.0;
                }
<<<<<<< HEAD
				
                /* ai -> aj interaction */
                if(raj < dr + sk_ai)
=======

                /* ai -> aj interaction */
                if (raj < dr + sk_ai)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    lij     = 1.0/(dr-sk_ai);
                    dlij    = 1.0;
                    raj_inv = 1.0/raj;
<<<<<<< HEAD
                    
                    if(raj>dr-sk_ai)
                    {
                        lij = raj_inv;
                        dlij = 0.0;
                    }
                    
                    lij2     = lij  * lij;
                    lij3     = lij2 * lij;
                    
                    uij      = 1.0/(dr+sk_ai);
                    uij2     = uij  * uij;
                    uij3     = uij2 * uij;
                    
                    diff2    = uij2-lij2;
                    
=======

                    if (raj > dr-sk_ai)
                    {
                        lij  = raj_inv;
                        dlij = 0.0;
                    }

                    lij2     = lij  * lij;
                    lij3     = lij2 * lij;

                    uij      = 1.0/(dr+sk_ai);
                    uij2     = uij  * uij;
                    uij3     = uij2 * uij;

                    diff2    = uij2-lij2;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    lij_inv  = gmx_invsqrt(lij2);
                    sk2      =  sk2_ai; /* sk2_ai = sk_ai * sk_ai in i loop above */
                    sk2_rinv = sk2*rinv;
                    prod     = 0.25 * sk2_rinv;
<<<<<<< HEAD
                    
                    /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                    log_term = log(uij*lij_inv);
                    
                    tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
                    
                    if(raj<sk_ai-dr)
                    {
                        tmp     = tmp + 2.0 * (raj_inv-lij);
                    }
                    
                    t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                    t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                    t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
                    
                    dadxj = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
                    
                    born->gpol_hct_work[k] += 0.5*tmp;
                }
                else 
=======

                    /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                    log_term = log(uij*lij_inv);

                    tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);

                    if (raj < sk_ai-dr)
                    {
                        tmp     = tmp + 2.0 * (raj_inv-lij);
                    }

                    t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                    t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                    t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                    dadxj = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/

                    born->gpol_hct_work[k] += 0.5*tmp;
                }
                else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    dadxj = 0.0;
                }
                fr->dadx[n++] = dadxi;
                fr->dadx[n++] = dadxj;
<<<<<<< HEAD
                
            }
        }
        
        /* Main part, no exclusions */
        for(j=nj1; j<nj2; j++)
        {       
            k = j%natoms;
            
=======

            }
        }

        /* Main part, no exclusions */
        for (j = nj1; j < nj2; j++)
        {
            k = j%natoms;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* load j atom coordinates */
            jx                = x[3*k];
            jy                = x[3*k+1];
            jz                = x[3*k+2];
<<<<<<< HEAD
            
            /* Calculate distance */
            dx                = ix - jx;      
            dy                = iy - jy;      
            dz                = iz - jz;      
            rsq               = dx*dx+dy*dy+dz*dz;
            
            /* Calculate 1/r and 1/r2 */
            rinv              = gmx_invsqrt(rsq);
            dr                = rsq*rinv;
            
            /* sk is precalculated in init_gb() */
            sk    = born->param[k];
            raj   = top->atomtypes.gb_radius[mdatoms->typeA[k]]-doffset; 
            
            /* aj -> ai interaction */
            if(rai < dr+sk)
            {
                lij       = 1.0/(dr-sk);
                dlij      = 1.0; 
                
                if(rai>dr-sk)
=======

            /* Calculate distance */
            dx                = ix - jx;
            dy                = iy - jy;
            dz                = iz - jz;
            rsq               = dx*dx+dy*dy+dz*dz;

            /* Calculate 1/r and 1/r2 */
            rinv              = gmx_invsqrt(rsq);
            dr                = rsq*rinv;

            /* sk is precalculated in init_gb() */
            sk    = born->param[k];
            raj   = top->atomtypes.gb_radius[mdatoms->typeA[k]]-doffset;

            /* aj -> ai interaction */
            if (rai < dr+sk)
            {
                lij       = 1.0/(dr-sk);
                dlij      = 1.0;

                if (rai > dr-sk)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    lij  = rai_inv;
                    dlij = 0.0;
                }
<<<<<<< HEAD
                
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                uij      = 1.0/(dr+sk);
                lij2     = lij  * lij;
                lij3     = lij2 * lij;
                uij2     = uij  * uij;
                uij3     = uij2 * uij;
<<<<<<< HEAD
                
                diff2    = uij2-lij2;
                
                lij_inv  = gmx_invsqrt(lij2);
                sk2      = sk*sk;
                sk2_rinv = sk2*rinv;	
                prod     = 0.25*sk2_rinv;
                
                log_term = log(uij*lij_inv);
                /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
                
                if(rai < sk-dr)
                {
                    tmp = tmp + 2.0 * (rai_inv-lij);
                }
                
                /* duij    = 1.0; */
                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr); 
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr); 
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv; 
                
                dadxi = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
                
                sum_ai += 0.5*tmp;
            }
            else 
=======

                diff2    = uij2-lij2;

                lij_inv  = gmx_invsqrt(lij2);
                sk2      = sk*sk;
                sk2_rinv = sk2*rinv;
                prod     = 0.25*sk2_rinv;

                log_term = log(uij*lij_inv);
                /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);

                if (rai < sk-dr)
                {
                    tmp = tmp + 2.0 * (rai_inv-lij);
                }

                /* duij    = 1.0; */
                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                dadxi = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/

                sum_ai += 0.5*tmp;
            }
            else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                dadxi = 0.0;
            }

            /* ai -> aj interaction */
<<<<<<< HEAD
            if(raj < dr + sk_ai)
=======
            if (raj < dr + sk_ai)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                lij     = 1.0/(dr-sk_ai);
                dlij    = 1.0;
                raj_inv = 1.0/raj;
<<<<<<< HEAD
                
                if(raj>dr-sk_ai)
                {
                    lij = raj_inv;
                    dlij = 0.0;
                }
                
                lij2     = lij  * lij;
                lij3     = lij2 * lij;
                
                uij      = 1.0/(dr+sk_ai);
                uij2     = uij  * uij;
                uij3     = uij2 * uij;
                
                diff2    = uij2-lij2;
                
=======

                if (raj > dr-sk_ai)
                {
                    lij  = raj_inv;
                    dlij = 0.0;
                }

                lij2     = lij  * lij;
                lij3     = lij2 * lij;

                uij      = 1.0/(dr+sk_ai);
                uij2     = uij  * uij;
                uij3     = uij2 * uij;

                diff2    = uij2-lij2;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                lij_inv  = gmx_invsqrt(lij2);
                sk2      =  sk2_ai; /* sk2_ai = sk_ai * sk_ai in i loop above */
                sk2_rinv = sk2*rinv;
                prod     = 0.25 * sk2_rinv;
<<<<<<< HEAD
                
                /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                log_term = log(uij*lij_inv);
                
                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
                
                if(raj<sk_ai-dr)
                {
                    tmp     = tmp + 2.0 * (raj_inv-lij);
                }
                
                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
                
                dadxj = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
                
                born->gpol_hct_work[k] += 0.5*tmp;
            }
            else 
=======

                /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                log_term = log(uij*lij_inv);

                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);

                if (raj < sk_ai-dr)
                {
                    tmp     = tmp + 2.0 * (raj_inv-lij);
                }

                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                dadxj = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/

                born->gpol_hct_work[k] += 0.5*tmp;
            }
            else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                dadxj = 0.0;
            }
            fr->dadx[n++] = dadxi;
<<<<<<< HEAD
            fr->dadx[n++] = dadxj;     
        }
        born->gpol_hct_work[i]+=sum_ai;
	}    
    
    /* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_hct_work, cr);
	}
	
    if(gb_algorithm==egbHCT)
    {
        /* HCT */
        for(i=0;i<natoms;i++)
        {
            if(born->use[i] != 0)
            {
                rai     = top->atomtypes.gb_radius[mdatoms->typeA[i]]-born->gb_doffset; 
                sum_ai  = 1.0/rai - born->gpol_hct_work[i];
                min_rad = rai + born->gb_doffset;
                rad     = 1.0/sum_ai; 
                
=======
            fr->dadx[n++] = dadxj;
        }
        born->gpol_hct_work[i] += sum_ai;
    }

    /* Parallel summations */
    if (PARTDECOMP(cr))
    {
        gmx_sum(natoms, born->gpol_hct_work, cr);
    }

    if (gb_algorithm == egbHCT)
    {
        /* HCT */
        for (i = 0; i < natoms; i++)
        {
            if (born->use[i] != 0)
            {
                rai     = top->atomtypes.gb_radius[mdatoms->typeA[i]]-born->gb_doffset;
                sum_ai  = 1.0/rai - born->gpol_hct_work[i];
                min_rad = rai + born->gb_doffset;
                rad     = 1.0/sum_ai;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                born->bRad[i]   = rad > min_rad ? rad : min_rad;
                fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
            }
        }
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        /* OBC */
        /* Calculate the radii */
<<<<<<< HEAD
        for(i=0;i<natoms;i++)
        {
            if(born->use[i] != 0)
            {
                rai        = top->atomtypes.gb_radius[mdatoms->typeA[i]];
                rai_inv2   = 1.0/rai;
                rai        = rai-doffset; 
=======
        for (i = 0; i < natoms; i++)
        {
            if (born->use[i] != 0)
            {
                rai        = top->atomtypes.gb_radius[mdatoms->typeA[i]];
                rai_inv2   = 1.0/rai;
                rai        = rai-doffset;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                rai_inv    = 1.0/rai;
                sum_ai     = rai * born->gpol_hct_work[i];
                sum_ai2    = sum_ai  * sum_ai;
                sum_ai3    = sum_ai2 * sum_ai;
<<<<<<< HEAD
                
                tsum    = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
                born->bRad[i] = rai_inv - tsum*rai_inv2;
                born->bRad[i] = 1.0 / born->bRad[i];
                
                fr->invsqrta[i]=gmx_invsqrt(born->bRad[i]);
                
                tchain  = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
                born->drobc[i] = (1.0-tsum*tsum)*tchain*rai_inv2;
            }
        }
    }   
=======

                tsum          = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
                born->bRad[i] = rai_inv - tsum*rai_inv2;
                born->bRad[i] = 1.0 / born->bRad[i];

                fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);

                tchain         = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
                born->drobc[i] = (1.0-tsum*tsum)*tchain*rai_inv2;
            }
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return 0;
}





int
genborn_allvsall_calc_chainrule(t_forcerec *           fr,
                                t_mdatoms *            mdatoms,
                                gmx_genborn_t *        born,
                                real *                 x,
                                real *                 f,
                                int                    gb_algorithm,
                                void *                 work)
{
<<<<<<< HEAD
	gmx_allvsallgb2_data_t *aadata;
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2;
	int        i,j,k,n;
    int        idx;
    int *      mask;
    
    real       ix,iy,iz;
    real       fix,fiy,fiz;
    real       jx,jy,jz;
    real       dx,dy,dz;
    real       tx,ty,tz;
    real       rbai,rbaj,fgb,fgb_ai,rbi;
    real *     rb;
    real *     dadx;
    
    natoms              = mdatoms->nr;
	ni0                 = mdatoms->start;
	ni1                 = mdatoms->start+mdatoms->homenr;
    dadx                = fr->dadx;
    
    aadata = (gmx_allvsallgb2_data_t *)work;

    n = 0;
    rb = born->work;
    
	/* Loop to get the proper form for the Born radius term */
	if(gb_algorithm==egbSTILL) 
	{
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = (2 * rbi * rbi * fr->dvda[i])/ONE_4PI_EPS0;
		}
	}
	else if(gb_algorithm==egbHCT) 
	{
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = rbi * rbi * fr->dvda[i];
		}
	}
	else if(gb_algorithm==egbOBC) 
	{
		for(idx=0;idx<natoms;idx++)
		{
			rbi   = born->bRad[idx];
			rb[idx] = rbi * rbi * born->drobc[idx] * fr->dvda[idx];
		}
	}

    for(i=ni0; i<ni1; i++)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */
		
		/* Load i atom data */
=======
    gmx_allvsallgb2_data_t *aadata;
    int                     natoms;
    int                     ni0, ni1;
    int                     nj0, nj1, nj2;
    int                     i, j, k, n;
    int                     idx;
    int              *      mask;

    real                    ix, iy, iz;
    real                    fix, fiy, fiz;
    real                    jx, jy, jz;
    real                    dx, dy, dz;
    real                    tx, ty, tz;
    real                    rbai, rbaj, fgb, fgb_ai, rbi;
    real              *     rb;
    real              *     dadx;

    natoms              = mdatoms->nr;
    ni0                 = mdatoms->start;
    ni1                 = mdatoms->start+mdatoms->homenr;
    dadx                = fr->dadx;

    aadata = (gmx_allvsallgb2_data_t *)work;

    n  = 0;
    rb = born->work;

    /* Loop to get the proper form for the Born radius term */
    if (gb_algorithm == egbSTILL)
    {
        for (i = 0; i < natoms; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = (2 * rbi * rbi * fr->dvda[i])/ONE_4PI_EPS0;
        }
    }
    else if (gb_algorithm == egbHCT)
    {
        for (i = 0; i < natoms; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = rbi * rbi * fr->dvda[i];
        }
    }
    else if (gb_algorithm == egbOBC)
    {
        for (idx = 0; idx < natoms; idx++)
        {
            rbi     = born->bRad[idx];
            rb[idx] = rbi * rbi * born->drobc[idx] * fr->dvda[idx];
        }
    }

    for (i = ni0; i < ni1; i++)
    {
        /* We assume shifts are NOT used for all-vs-all interactions */

        /* Load i atom data */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        ix                = x[3*i];
        iy                = x[3*i+1];
        iz                = x[3*i+2];

        fix               = 0;
        fiy               = 0;
        fiz               = 0;
<<<<<<< HEAD
        
        rbai              = rb[i];
        
		/* Load limits for loop over neighbors */
		nj0              = aadata->jindex_gb[3*i];
		nj1              = aadata->jindex_gb[3*i+1];
		nj2              = aadata->jindex_gb[3*i+2];
        
        mask             = aadata->exclusion_mask_gb[i];
                
        /* Prologue part, including exclusion mask */
        for(j=nj0; j<nj1; j++,mask++)
        {          
            if(*mask!=0)
            {
                k = j%natoms;
                
=======

        rbai              = rb[i];

        /* Load limits for loop over neighbors */
        nj0              = aadata->jindex_gb[3*i];
        nj1              = aadata->jindex_gb[3*i+1];
        nj2              = aadata->jindex_gb[3*i+2];

        mask             = aadata->exclusion_mask_gb[i];

        /* Prologue part, including exclusion mask */
        for (j = nj0; j < nj1; j++, mask++)
        {
            if (*mask != 0)
            {
                k = j%natoms;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                /* load j atom coordinates */
                jx                = x[3*k];
                jy                = x[3*k+1];
                jz                = x[3*k+2];
<<<<<<< HEAD
                
                /* Calculate distance */
                dx                = ix - jx;      
                dy                = iy - jy;      
                dz                = iz - jz;   
                
                rbaj              = rb[k];
                
                fgb     = rbai*dadx[n++]; 
                fgb_ai  = rbaj*dadx[n++];
                
                /* Total force between ai and aj is the sum of ai->aj and aj->ai */
                fgb     = fgb + fgb_ai;
                
                tx      = fgb * dx;
                ty      = fgb * dy;
                tz      = fgb * dz;
                
                fix     = fix + tx;
                fiy     = fiy + ty;
                fiz     = fiz + tz;
                
=======

                /* Calculate distance */
                dx                = ix - jx;
                dy                = iy - jy;
                dz                = iz - jz;

                rbaj              = rb[k];

                fgb     = rbai*dadx[n++];
                fgb_ai  = rbaj*dadx[n++];

                /* Total force between ai and aj is the sum of ai->aj and aj->ai */
                fgb     = fgb + fgb_ai;

                tx      = fgb * dx;
                ty      = fgb * dy;
                tz      = fgb * dz;

                fix     = fix + tx;
                fiy     = fiy + ty;
                fiz     = fiz + tz;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                /* Update force on atom aj */
                f[3*k]   = f[3*k] - tx;
                f[3*k+1] = f[3*k+1] - ty;
                f[3*k+2] = f[3*k+2] - tz;
            }
        }
<<<<<<< HEAD
        
        /* Main part, no exclusions */
        for(j=nj1; j<nj2; j++)
        {       
            k = j%natoms;
            
=======

        /* Main part, no exclusions */
        for (j = nj1; j < nj2; j++)
        {
            k = j%natoms;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* load j atom coordinates */
            jx                = x[3*k];
            jy                = x[3*k+1];
            jz                = x[3*k+2];
<<<<<<< HEAD
            
            /* Calculate distance */
            dx                = ix - jx;      
            dy                = iy - jy;      
            dz                = iz - jz;   
            
            rbaj              = rb[k];
            
            fgb     = rbai*dadx[n++]; 
=======

            /* Calculate distance */
            dx                = ix - jx;
            dy                = iy - jy;
            dz                = iz - jz;

            rbaj              = rb[k];

            fgb     = rbai*dadx[n++];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            fgb_ai  = rbaj*dadx[n++];

            /* Total force between ai and aj is the sum of ai->aj and aj->ai */
            fgb     = fgb + fgb_ai;
<<<<<<< HEAD
            
            tx      = fgb * dx;
            ty      = fgb * dy;
            tz      = fgb * dz;
            
            fix     = fix + tx;
            fiy     = fiy + ty;
            fiz     = fiz + tz;
            
=======

            tx      = fgb * dx;
            ty      = fgb * dy;
            tz      = fgb * dz;

            fix     = fix + tx;
            fiy     = fiy + ty;
            fiz     = fiz + tz;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* Update force on atom aj */
            f[3*k]   = f[3*k] - tx;
            f[3*k+1] = f[3*k+1] - ty;
            f[3*k+2] = f[3*k+2] - tz;
        }
<<<<<<< HEAD
		/* Update force and shift forces on atom ai */
		f[3*i]   = f[3*i] + fix;
		f[3*i+1] = f[3*i+1] + fiy;
		f[3*i+2] = f[3*i+2] + fiz;
	}    
        
	return 0;
}



=======
        /* Update force and shift forces on atom ai */
        f[3*i]   = f[3*i] + fix;
        f[3*i+1] = f[3*i+1] + fiy;
        f[3*i+2] = f[3*i+2] + fiz;
    }

    return 0;
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
