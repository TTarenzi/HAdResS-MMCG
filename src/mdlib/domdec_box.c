<<<<<<< HEAD
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2008
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

#include "typedefs.h"
#include "vec.h"
#include "pbc.h"
#include "domdec.h"
#include "domdec_network.h"
#include "nsgrid.h"
#include "network.h"

<<<<<<< HEAD
static void calc_cgcm_av_stddev(t_block *cgs,int n,rvec *x,rvec av,rvec stddev,
                                t_commrec *cr_sum)
{
    int  *cgindex;
    dvec s1,s2;
    double buf[7];
    int  cg,d,k0,k1,k,nrcg;
    real inv_ncg;
    rvec cg_cm;
=======
static void calc_cgcm_av_stddev(t_block *cgs, int n, rvec *x, rvec av, rvec stddev,
                                t_commrec *cr_sum)
{
    int   *cgindex;
    dvec   s1, s2;
    double buf[7];
    int    cg, d, k0, k1, k, nrcg;
    real   inv_ncg;
    rvec   cg_cm;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    clear_dvec(s1);
    clear_dvec(s2);

    cgindex = cgs->index;
<<<<<<< HEAD
    for(cg=0; cg<n; cg++)
=======
    for (cg = 0; cg < n; cg++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        k0      = cgindex[cg];
        k1      = cgindex[cg+1];
        nrcg    = k1 - k0;
        if (nrcg == 1)
        {
<<<<<<< HEAD
            copy_rvec(x[k0],cg_cm);
=======
            copy_rvec(x[k0], cg_cm);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        else
        {
            inv_ncg = 1.0/nrcg;
<<<<<<< HEAD
            
            clear_rvec(cg_cm);
            for(k=k0; (k<k1); k++)
            {
                rvec_inc(cg_cm,x[k]);
            }
            for(d=0; (d<DIM); d++)
=======

            clear_rvec(cg_cm);
            for (k = k0; (k < k1); k++)
            {
                rvec_inc(cg_cm, x[k]);
            }
            for (d = 0; (d < DIM); d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                cg_cm[d] *= inv_ncg;
            }
        }
<<<<<<< HEAD
        for(d=0; d<DIM; d++)
=======
        for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            s1[d] += cg_cm[d];
            s2[d] += cg_cm[d]*cg_cm[d];
        }
    }

    if (cr_sum != NULL)
    {
<<<<<<< HEAD
        for(d=0; d<DIM; d++)
=======
        for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            buf[d]     = s1[d];
            buf[DIM+d] = s2[d];
        }
        buf[6] = n;
<<<<<<< HEAD
        gmx_sumd(7,buf,cr_sum);
        for(d=0; d<DIM; d++)
=======
        gmx_sumd(7, buf, cr_sum);
        for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            s1[d] = buf[d];
            s2[d] = buf[DIM+d];
        }
        n = (int)(buf[6] + 0.5);
    }

<<<<<<< HEAD
    dsvmul(1.0/n,s1,s1);
    dsvmul(1.0/n,s2,s2);

    for(d=0; d<DIM; d++)
=======
    dsvmul(1.0/n, s1, s1);
    dsvmul(1.0/n, s2, s2);

    for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        av[d]     = s1[d];
        stddev[d] = sqrt(s2[d] - s1[d]*s1[d]);
    }
}

<<<<<<< HEAD
static void set_tric_dir(ivec *dd_nc,gmx_ddbox_t *ddbox,matrix box)
{
    int  npbcdim,d,i,j;
    rvec *v,*normal;
    real dep,inv_skew_fac2;
    
    npbcdim = ddbox->npbcdim;
    normal  = ddbox->normal;
    for(d=0; d<DIM; d++)
    {
        ddbox->tric_dir[d] = 0;
        for(j=d+1; j<npbcdim; j++)
=======
static void set_tric_dir(ivec *dd_nc, gmx_ddbox_t *ddbox, matrix box)
{
    int   npbcdim, d, i, j;
    rvec *v, *normal;
    real  dep, inv_skew_fac2;

    npbcdim = ddbox->npbcdim;
    normal  = ddbox->normal;
    for (d = 0; d < DIM; d++)
    {
        ddbox->tric_dir[d] = 0;
        for (j = d+1; j < npbcdim; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            if (box[j][d] != 0)
            {
                ddbox->tric_dir[d] = 1;
                if (dd_nc != NULL && (*dd_nc)[j] > 1 && (*dd_nc)[d] == 1)
                {
<<<<<<< HEAD
                    gmx_fatal(FARGS,"Domain decomposition has not been implemented for box vectors that have non-zero components in directions that do not use domain decomposition: ncells = %d %d %d, box vector[%d] = %f %f %f",
                              dd_nc[XX],dd_nc[YY],dd_nc[ZZ],
                              j+1,box[j][XX],box[j][YY],box[j][ZZ]);
                }
            }
        }
        
=======
                    gmx_fatal(FARGS, "Domain decomposition has not been implemented for box vectors that have non-zero components in directions that do not use domain decomposition: ncells = %d %d %d, box vector[%d] = %f %f %f",
                              dd_nc[XX], dd_nc[YY], dd_nc[ZZ],
                              j+1, box[j][XX], box[j][YY], box[j][ZZ]);
                }
            }
        }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Convert box vectors to orthogonal vectors for this dimension,
         * for use in distance calculations.
         * Set the trilinic skewing factor that translates
         * the thickness of a slab perpendicular to this dimension
         * into the real thickness of the slab.
         */
        if (ddbox->tric_dir[d])
        {
            inv_skew_fac2 = 1;
<<<<<<< HEAD
            v = ddbox->v[d];
            if (d == XX || d == YY)
            {
                /* Normalize such that the "diagonal" is 1 */
                svmul(1/box[d+1][d+1],box[d+1],v[d+1]);
                for(i=0; i<d; i++)
=======
            v             = ddbox->v[d];
            if (d == XX || d == YY)
            {
                /* Normalize such that the "diagonal" is 1 */
                svmul(1/box[d+1][d+1], box[d+1], v[d+1]);
                for (i = 0; i < d; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    v[d+1][i] = 0;
                }
                inv_skew_fac2 += sqr(v[d+1][d]);
                if (d == XX)
                {
                    /* Normalize such that the "diagonal" is 1 */
<<<<<<< HEAD
                    svmul(1/box[d+2][d+2],box[d+2],v[d+2]);
                    for(i=0; i<d; i++)
=======
                    svmul(1/box[d+2][d+2], box[d+2], v[d+2]);
                    for (i = 0; i < d; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        v[d+2][i] = 0;
                    }
                    /* Make vector [d+2] perpendicular to vector [d+1],
                     * this does not affect the normalization.
                     */
<<<<<<< HEAD
                    dep = iprod(v[d+1],v[d+2])/norm2(v[d+1]);
                    for(i=0; i<DIM; i++)
=======
                    dep = iprod(v[d+1], v[d+2])/norm2(v[d+1]);
                    for (i = 0; i < DIM; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        v[d+2][i] -= dep*v[d+1][i];
                    }
                    inv_skew_fac2 += sqr(v[d+2][d]);
<<<<<<< HEAD
                    
                    cprod(v[d+1],v[d+2],normal[d]);
=======

                    cprod(v[d+1], v[d+2], normal[d]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
                else
                {
                    /* cross product with (1,0,0) */
                    normal[d][XX] =  0;
                    normal[d][YY] =  v[d+1][ZZ];
                    normal[d][ZZ] = -v[d+1][YY];
                }
                if (debug)
                {
<<<<<<< HEAD
                    fprintf(debug,"box[%d]  %.3f %.3f %.3f\n",
                            d,box[d][XX],box[d][YY],box[d][ZZ]);
                    for(i=d+1; i<DIM; i++)
                    {
                        fprintf(debug,"  v[%d]  %.3f %.3f %.3f\n",
                                i,v[i][XX],v[i][YY],v[i][ZZ]);
=======
                    fprintf(debug, "box[%d]  %.3f %.3f %.3f\n",
                            d, box[d][XX], box[d][YY], box[d][ZZ]);
                    for (i = d+1; i < DIM; i++)
                    {
                        fprintf(debug, "  v[%d]  %.3f %.3f %.3f\n",
                                i, v[i][XX], v[i][YY], v[i][ZZ]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    }
                }
            }
            ddbox->skew_fac[d] = 1.0/sqrt(inv_skew_fac2);
            /* Set the normal vector length to skew_fac */
            dep = ddbox->skew_fac[d]/norm(normal[d]);
<<<<<<< HEAD
            svmul(dep,normal[d],normal[d]);

            if (debug)
            {
                fprintf(debug,"skew_fac[%d] = %f\n",d,ddbox->skew_fac[d]);
                fprintf(debug,"normal[%d]  %.3f %.3f %.3f\n",
                        d,normal[d][XX],normal[d][YY],normal[d][ZZ]);
=======
            svmul(dep, normal[d], normal[d]);

            if (debug)
            {
                fprintf(debug, "skew_fac[%d] = %f\n", d, ddbox->skew_fac[d]);
                fprintf(debug, "normal[%d]  %.3f %.3f %.3f\n",
                        d, normal[d][XX], normal[d][YY], normal[d][ZZ]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
        else
        {
            ddbox->skew_fac[d] = 1;
<<<<<<< HEAD
            
            for(i=0; i<DIM; i++)
            {
                clear_rvec(ddbox->v[d][i]);
                ddbox->v[d][i][i] = 1;
            }   
=======

            for (i = 0; i < DIM; i++)
            {
                clear_rvec(ddbox->v[d][i]);
                ddbox->v[d][i][i] = 1;
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            clear_rvec(normal[d]);
            normal[d][d] = 1;
        }
    }
}

<<<<<<< HEAD
static void low_set_ddbox(t_inputrec *ir,ivec *dd_nc,matrix box,
                          gmx_bool bCalcUnboundedSize,int ncg,t_block *cgs,rvec *x,
                          t_commrec *cr_sum,
                          gmx_ddbox_t *ddbox)
{
    rvec av,stddev;
    real b0,b1;
=======
static void low_set_ddbox(t_inputrec *ir, ivec *dd_nc, matrix box,
                          gmx_bool bCalcUnboundedSize, int ncg, t_block *cgs, rvec *x,
                          t_commrec *cr_sum,
                          gmx_ddbox_t *ddbox)
{
    rvec av, stddev;
    real b0, b1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    int  d;

    ddbox->npbcdim     = ePBC2npbcdim(ir->ePBC);
    ddbox->nboundeddim = inputrec2nboundeddim(ir);

<<<<<<< HEAD
    for(d=0; d<ddbox->nboundeddim; d++)
=======
    for (d = 0; d < ddbox->nboundeddim; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        ddbox->box0[d]     = 0;
        ddbox->box_size[d] = box[d][d];
    }

    if (ddbox->nboundeddim < DIM && bCalcUnboundedSize)
    {
<<<<<<< HEAD
        calc_cgcm_av_stddev(cgs,ncg,x,av,stddev,cr_sum);
=======
        calc_cgcm_av_stddev(cgs, ncg, x, av, stddev, cr_sum);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* GRID_STDDEV_FAC * stddev
         * gives a uniform load for a rectangular block of cg's.
         * For a sphere it is not a bad approximation for 4x1x1 up to 4x2x2.
         */
<<<<<<< HEAD
        for(d=ddbox->nboundeddim; d<DIM; d++)
=======
        for (d = ddbox->nboundeddim; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            b0 = av[d] - GRID_STDDEV_FAC*stddev[d];
            b1 = av[d] + GRID_STDDEV_FAC*stddev[d];
            if (debug)
            {
<<<<<<< HEAD
                fprintf(debug,"Setting global DD grid boundaries to %f - %f\n",
                        b0,b1);
=======
                fprintf(debug, "Setting global DD grid boundaries to %f - %f\n",
                        b0, b1);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            ddbox->box0[d]     = b0;
            ddbox->box_size[d] = b1 - b0;
        }
    }

<<<<<<< HEAD
    set_tric_dir(dd_nc,ddbox,box);
}

void set_ddbox(gmx_domdec_t *dd,gmx_bool bMasterState,t_commrec *cr_sum,
               t_inputrec *ir,matrix box,
               gmx_bool bCalcUnboundedSize,t_block *cgs,rvec *x,
=======
    set_tric_dir(dd_nc, ddbox, box);
}

void set_ddbox(gmx_domdec_t *dd, gmx_bool bMasterState, t_commrec *cr_sum,
               t_inputrec *ir, matrix box,
               gmx_bool bCalcUnboundedSize, t_block *cgs, rvec *x,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
               gmx_ddbox_t *ddbox)
{
    if (!bMasterState || DDMASTER(dd))
    {
<<<<<<< HEAD
        low_set_ddbox(ir,&dd->nc,box,bCalcUnboundedSize,
                      bMasterState ? cgs->nr : dd->ncg_home,cgs,x,
=======
        low_set_ddbox(ir, &dd->nc, box, bCalcUnboundedSize,
                      bMasterState ? cgs->nr : dd->ncg_home, cgs, x,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                      bMasterState ? NULL : cr_sum,
                      ddbox);
    }

    if (bMasterState)
    {
<<<<<<< HEAD
        dd_bcast(dd,sizeof(gmx_ddbox_t),ddbox);
    }
}

void set_ddbox_cr(t_commrec *cr,ivec *dd_nc,
                  t_inputrec *ir,matrix box,t_block *cgs,rvec *x,
=======
        dd_bcast(dd, sizeof(gmx_ddbox_t), ddbox);
    }
}

void set_ddbox_cr(t_commrec *cr, ivec *dd_nc,
                  t_inputrec *ir, matrix box, t_block *cgs, rvec *x,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                  gmx_ddbox_t *ddbox)
{
    if (MASTER(cr))
    {
<<<<<<< HEAD
        low_set_ddbox(ir,dd_nc,box,TRUE,cgs->nr,cgs,x,NULL,ddbox);
    }

    gmx_bcast(sizeof(gmx_ddbox_t),ddbox,cr);
=======
        low_set_ddbox(ir, dd_nc, box, TRUE, cgs->nr, cgs, x, NULL, ddbox);
    }

    gmx_bcast(sizeof(gmx_ddbox_t), ddbox, cr);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}
