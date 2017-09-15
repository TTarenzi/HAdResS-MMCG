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

#include <stdio.h>
<<<<<<< HEAD
=======
#include <assert.h>
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "domdec.h"
#include "network.h"
#include "perf_est.h"
#include "physics.h"
#include "smalloc.h"
#include "typedefs.h"
#include "vec.h"
#include "names.h"

/* Margin for setting up the DD grid */
#define DD_GRID_MARGIN_PRES_SCALE 1.05

<<<<<<< HEAD
static int factorize(int n,int **fac,int **mfac)
{
    int d,ndiv;
=======
static int factorize(int n, int **fac, int **mfac)
{
    int d, ndiv;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (n <= 0)
    {
        gmx_fatal(FARGS, "Can only factorize positive integers.");
    }

    /* Decompose n in factors */
<<<<<<< HEAD
    snew(*fac,n/2);
    snew(*mfac,n/2);
    d = 2;
=======
    snew(*fac, n/2);
    snew(*mfac, n/2);
    d    = 2;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    ndiv = 0;
    while (n > 1)
    {
        while (n % d == 0)
        {
            if (ndiv == 0 || (*fac)[ndiv-1] != d)
            {
                ndiv++;
                (*fac)[ndiv-1] = d;
            }
            (*mfac)[ndiv-1]++;
            n /= d;
        }
        d++;
    }
<<<<<<< HEAD
	
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return ndiv;
}

static gmx_bool largest_divisor(int n)
{
<<<<<<< HEAD
    int ndiv,*div,*mdiv,ldiv;

    ndiv = factorize(n,&div,&mdiv);
=======
    int ndiv, *div, *mdiv, ldiv;

    ndiv = factorize(n, &div, &mdiv);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    ldiv = div[ndiv-1];
    sfree(div);
    sfree(mdiv);

    return ldiv;
}

<<<<<<< HEAD
static int lcd(int n1,int n2)
{
    int d,i;
    
    d = 1;
    for(i=2; (i<=n1 && i<=n2); i++)
=======
static int lcd(int n1, int n2)
{
    int d, i;

    d = 1;
    for (i = 2; (i <= n1 && i <= n2); i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (n1 % i == 0 && n2 % i == 0)
        {
            d = i;
        }
    }
<<<<<<< HEAD
    
  return d;
}

static gmx_bool fits_pme_ratio(int nnodes,int npme,float ratio)
{
    return ((double)npme/(double)nnodes > 0.95*ratio); 
}

static gmx_bool fits_pp_pme_perf(FILE *fplog,
                             t_inputrec *ir,matrix box,gmx_mtop_t *mtop,
                             int nnodes,int npme,float ratio)
{
    int ndiv,*div,*mdiv,ldiv;
    int npp_root3,npme_root2;

    ndiv = factorize(nnodes-npme,&div,&mdiv);
=======

    return d;
}

static gmx_bool fits_pme_ratio(int nnodes, int npme, float ratio)
{
    return ((double)npme/(double)nnodes > 0.95*ratio);
}

static gmx_bool fits_pp_pme_perf(FILE *fplog,
                                 t_inputrec *ir, matrix box, gmx_mtop_t *mtop,
                                 int nnodes, int npme, float ratio)
{
    int ndiv, *div, *mdiv, ldiv;
    int npp_root3, npme_root2;

    ndiv = factorize(nnodes-npme, &div, &mdiv);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    ldiv = div[ndiv-1];
    sfree(div);
    sfree(mdiv);

<<<<<<< HEAD
    npp_root3  = (int)(pow(nnodes-npme,1.0/3.0) + 0.5);
=======
    npp_root3  = (int)(pow(nnodes-npme, 1.0/3.0) + 0.5);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    npme_root2 = (int)(sqrt(npme) + 0.5);

    /* The check below gives a reasonable division:
     * factor 5 allowed at 5 or more PP nodes,
     * factor 7 allowed at 49 or more PP nodes.
     */
    if (ldiv > 3 + npp_root3)
    {
        return FALSE;
    }

    /* Check if the number of PP and PME nodes have a reasonable sized
     * denominator in common, such that we can use 2D PME decomposition
     * when required (which requires nx_pp == nx_pme).
     * The factor of 2 allows for a maximum ratio of 2^2=4
     * between nx_pme and ny_pme.
     */
<<<<<<< HEAD
    if (lcd(nnodes-npme,npme)*2 < npme_root2)
=======
    if (lcd(nnodes-npme, npme)*2 < npme_root2)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return FALSE;
    }

    /* Does this division gives a reasonable PME load? */
<<<<<<< HEAD
    return fits_pme_ratio(nnodes,npme,ratio);
}

static int guess_npme(FILE *fplog,gmx_mtop_t *mtop,t_inputrec *ir,matrix box,
					  int nnodes)
{
	float ratio;
	int  npme,nkx,nky;
	t_inputrec ir_try;
	
	ratio = pme_load_estimate(mtop,ir,box);
	
	if (fplog)
    {
		fprintf(fplog,"Guess for relative PME load: %.2f\n",ratio);
    }
	
	/* We assume the optimal node ratio is close to the load ratio.
	 * The communication load is neglected,
	 * but (hopefully) this will balance out between PP and PME.
	 */
	
    if (!fits_pme_ratio(nnodes,nnodes/2,ratio))
=======
    return fits_pme_ratio(nnodes, npme, ratio);
}

static int guess_npme(FILE *fplog, gmx_mtop_t *mtop, t_inputrec *ir, matrix box,
                      int nnodes)
{
    float      ratio;
    int        npme, nkx, nky;
    t_inputrec ir_try;

    ratio = pme_load_estimate(mtop, ir, box);

    if (fplog)
    {
        fprintf(fplog, "Guess for relative PME load: %.2f\n", ratio);
    }

    /* We assume the optimal node ratio is close to the load ratio.
     * The communication load is neglected,
     * but (hopefully) this will balance out between PP and PME.
     */

    if (!fits_pme_ratio(nnodes, nnodes/2, ratio))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* We would need more than nnodes/2 PME only nodes,
         * which is not possible. Since the PME load is very high,
         * we will not loose much performance when all nodes do PME.
         */

        return 0;
    }

    /* First try to find npme as a factor of nnodes up to nnodes/3.
     * We start with a minimum PME node fraction of 1/16
     * and avoid ratios which lead to large prime factors in nnodes-npme.
     */
    npme = (nnodes + 15)/16;
<<<<<<< HEAD
    while (npme <= nnodes/3) {
=======
    while (npme <= nnodes/3)
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (nnodes % npme == 0)
        {
            /* Note that fits_perf might change the PME grid,
             * in the current implementation it does not.
             */
<<<<<<< HEAD
            if (fits_pp_pme_perf(fplog,ir,box,mtop,nnodes,npme,ratio))
			{
				break;
			}
=======
            if (fits_pp_pme_perf(fplog, ir, box, mtop, nnodes, npme, ratio))
            {
                break;
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        npme++;
    }
    if (npme > nnodes/3)
    {
        /* Try any possible number for npme */
        npme = 1;
        while (npme <= nnodes/2)
        {
            /* Note that fits_perf may change the PME grid */
<<<<<<< HEAD
            if (fits_pp_pme_perf(fplog,ir,box,mtop,nnodes,npme,ratio))
=======
            if (fits_pp_pme_perf(fplog, ir, box, mtop, nnodes, npme, ratio))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                break;
            }
            npme++;
        }
    }
    if (npme > nnodes/2)
    {
<<<<<<< HEAD
        gmx_fatal(FARGS,"Could not find an appropriate number of separate PME nodes. i.e. >= %5f*#nodes (%d) and <= #nodes/2 (%d) and reasonable performance wise (grid_x=%d, grid_y=%d).\n"
                  "Use the -npme option of mdrun or change the number of processors or the PME grid dimensions, see the manual for details.",
                  ratio,(int)(0.95*ratio*nnodes+0.5),nnodes/2,ir->nkx,ir->nky);
=======
        gmx_fatal(FARGS, "Could not find an appropriate number of separate PME nodes. i.e. >= %5f*#nodes (%d) and <= #nodes/2 (%d) and reasonable performance wise (grid_x=%d, grid_y=%d).\n"
                  "Use the -npme option of mdrun or change the number of processors or the PME grid dimensions, see the manual for details.",
                  ratio, (int)(0.95*ratio*nnodes+0.5), nnodes/2, ir->nkx, ir->nky);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Keep the compiler happy */
        npme = 0;
    }
    else
    {
        if (fplog)
        {
            fprintf(fplog,
                    "Will use %d particle-particle and %d PME only nodes\n"
                    "This is a guess, check the performance at the end of the log file\n",
<<<<<<< HEAD
                    nnodes-npme,npme);
        }
        fprintf(stderr,"\n"
                "Will use %d particle-particle and %d PME only nodes\n"
                "This is a guess, check the performance at the end of the log file\n",
                nnodes-npme,npme);
    }
    
    return npme;
}

static int div_up(int n,int f)
=======
                    nnodes-npme, npme);
        }
        fprintf(stderr, "\n"
                "Will use %d particle-particle and %d PME only nodes\n"
                "This is a guess, check the performance at the end of the log file\n",
                nnodes-npme, npme);
    }

    return npme;
}

static int div_up(int n, int f)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    return (n + f - 1)/f;
}

<<<<<<< HEAD
real comm_box_frac(ivec dd_nc,real cutoff,gmx_ddbox_t *ddbox)
{
    int  i,j,k,npp;
    rvec bt,nw;
    real comm_vol;

    for(i=0; i<DIM; i++)
=======
real comm_box_frac(ivec dd_nc, real cutoff, gmx_ddbox_t *ddbox)
{
    int  i, j, k, npp;
    rvec bt, nw;
    real comm_vol;

    for (i = 0; i < DIM; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        bt[i] = ddbox->box_size[i]*ddbox->skew_fac[i];
        nw[i] = dd_nc[i]*cutoff/bt[i];
    }

<<<<<<< HEAD
    npp = 1;
    comm_vol = 0;
    for(i=0; i<DIM; i++)
    {
        if (dd_nc[i] > 1)
        {
            npp *= dd_nc[i];
            comm_vol += nw[i];
            for(j=i+1; j<DIM; j++)
=======
    npp      = 1;
    comm_vol = 0;
    for (i = 0; i < DIM; i++)
    {
        if (dd_nc[i] > 1)
        {
            npp      *= dd_nc[i];
            comm_vol += nw[i];
            for (j = i+1; j < DIM; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                if (dd_nc[j] > 1)
                {
                    comm_vol += nw[i]*nw[j]*M_PI/4;
<<<<<<< HEAD
                    for(k=j+1; k<DIM; k++)
=======
                    for (k = j+1; k < DIM; k++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        if (dd_nc[k] > 1)
                        {
                            comm_vol += nw[i]*nw[j]*nw[k]*M_PI/6;
                        }
                    }
                }
            }
        }
    }
<<<<<<< HEAD
   
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return comm_vol;
}

static gmx_bool inhomogeneous_z(const t_inputrec *ir)
{
<<<<<<< HEAD
    return ((EEL_PME(ir->coulombtype) || ir->coulombtype==eelEWALD) &&
            ir->ePBC==epbcXYZ && ir->ewald_geometry==eewg3DC);
=======
    return ((EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD) &&
            ir->ePBC == epbcXYZ && ir->ewald_geometry == eewg3DC);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

/* Avoid integer overflows */
static float comm_pme_cost_vol(int npme, int a, int b, int c)
{
    float comm_vol;

<<<<<<< HEAD
    comm_vol = npme - 1;
=======
    comm_vol  = npme - 1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    comm_vol *= npme;
    comm_vol *= div_up(a, npme);
    comm_vol *= div_up(b, npme);
    comm_vol *= c;
    return comm_vol;
}

<<<<<<< HEAD
static float comm_cost_est(gmx_domdec_t *dd,real limit,real cutoff,
                           matrix box,gmx_ddbox_t *ddbox,
                           int natoms,t_inputrec *ir,
                           float pbcdxr,
                           int npme_tot,ivec nc)
{
    ivec npme={1,1,1};
    int  i,j,k,nk,overlap;
    rvec bt;
    float comm_vol,comm_vol_xf,comm_pme,cost_pbcdx;
=======
static float comm_cost_est(gmx_domdec_t *dd, real limit, real cutoff,
                           matrix box, gmx_ddbox_t *ddbox,
                           int natoms, t_inputrec *ir,
                           float pbcdxr,
                           int npme_tot, ivec nc)
{
    ivec  npme = {1, 1, 1};
    int   i, j, k, nk, overlap;
    rvec  bt;
    float comm_vol, comm_vol_xf, comm_pme, cost_pbcdx;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* This is the cost of a pbc_dx call relative to the cost
     * of communicating the coordinate and force of an atom.
     * This will be machine dependent.
     * These factors are for x86 with SMP or Infiniband.
     */
    float pbcdx_rect_fac = 0.1;
    float pbcdx_tric_fac = 0.2;
    float temp;
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Check the DD algorithm restrictions */
    if ((ir->ePBC == epbcXY && ir->nwall < 2 && nc[ZZ] > 1) ||
        (ir->ePBC == epbcSCREW && (nc[XX] == 1 || nc[YY] > 1 || nc[ZZ] > 1)))
    {
        return -1;
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (inhomogeneous_z(ir) && nc[ZZ] > 1)
    {
        return -1;
    }

<<<<<<< HEAD
    /* Check if the triclinic requirements are met */
    for(i=0; i<DIM; i++)
    {
        for(j=i+1; j<ddbox->npbcdim; j++)
=======
    assert(ddbox->npbcdim <= DIM);

    /* Check if the triclinic requirements are met */
    for (i = 0; i < DIM; i++)
    {
        for (j = i+1; j < ddbox->npbcdim; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            if (box[j][i] != 0 || ir->deform[j][i] != 0 ||
                (ir->epc != epcNO && ir->compress[j][i] != 0))
            {
                if (nc[j] > 1 && nc[i] == 1)
                {
                    return -1;
                }
            }
        }
    }
<<<<<<< HEAD
    
    for(i=0; i<DIM; i++)
    {
        bt[i] = ddbox->box_size[i]*ddbox->skew_fac[i];
        
=======

    for (i = 0; i < DIM; i++)
    {
        bt[i] = ddbox->box_size[i]*ddbox->skew_fac[i];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Without PBC there are no cell size limits with 2 cells */
        if (!(i >= ddbox->npbcdim && nc[i] <= 2) && bt[i] < nc[i]*limit)
        {
            return -1;
        }
    }

    if (npme_tot > 1)
    {
        /* The following choices should match those
         * in init_domain_decomposition in domdec.c.
         */
        if (nc[XX] == 1 && nc[YY] > 1)
        {
            npme[XX] = 1;
            npme[YY] = npme_tot;
        }
        else if (nc[YY] == 1)
        {
            npme[XX] = npme_tot;
            npme[YY] = 1;
        }
        else
        {
            /* Will we use 1D or 2D PME decomposition? */
            npme[XX] = (npme_tot % nc[XX] == 0) ? nc[XX] : npme_tot;
            npme[YY] = npme_tot/npme[XX];
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* When two dimensions are (nearly) equal, use more cells
     * for the smallest index, so the decomposition does not
     * depend sensitively on the rounding of the box elements.
     */
<<<<<<< HEAD
    for(i=0; i<DIM; i++)
    {
        for(j=i+1; j<DIM; j++)
=======
    for (i = 0; i < DIM; i++)
    {
        for (j = i+1; j < DIM; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* Check if the box size is nearly identical,
             * in that case we prefer nx > ny  and ny > nz.
             */
            if (fabs(bt[j] - bt[i]) < 0.01*bt[i] && nc[j] > nc[i])
            {
                /* The XX/YY check is a bit compact. If nc[YY]==npme[YY]
<<<<<<< HEAD
             * this means the swapped nc has nc[XX]==npme[XX],
             * and we can also swap X and Y for PME.
             */
=======
                 * this means the swapped nc has nc[XX]==npme[XX],
                 * and we can also swap X and Y for PME.
                 */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                /* Check if dimension i and j are equivalent for PME.
                 * For x/y: if nc[YY]!=npme[YY], we can not swap x/y
                 * For y/z: we can not have PME decomposition in z
                 */
                if (npme_tot <= 1 ||
                    !((i == XX && j == YY && nc[YY] != npme[YY]) ||
                      (i == YY && j == ZZ && npme[YY] > 1)))
                {
                    return -1;
                }
            }
        }
    }

    /* This function determines only half of the communication cost.
     * All PP, PME and PP-PME communication is symmetric
     * and the "back"-communication cost is identical to the forward cost.
     */
<<<<<<< HEAD
    
    comm_vol = comm_box_frac(nc,cutoff,ddbox);

    comm_pme = 0;
    for(i=0; i<2; i++)
=======

    comm_vol = comm_box_frac(nc, cutoff, ddbox);

    comm_pme = 0;
    for (i = 0; i < 2; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* Determine the largest volume for PME x/f redistribution */
        if (nc[i] % npme[i] != 0)
        {
            if (nc[i] > npme[i])
            {
<<<<<<< HEAD
                comm_vol_xf = (npme[i]==2 ? 1.0/3.0 : 0.5);
            }
            else
            {
                comm_vol_xf = 1.0 - lcd(nc[i],npme[i])/(double)npme[i];
=======
                comm_vol_xf = (npme[i] == 2 ? 1.0/3.0 : 0.5);
            }
            else
            {
                comm_vol_xf = 1.0 - lcd(nc[i], npme[i])/(double)npme[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            comm_pme += 3*natoms*comm_vol_xf;
        }

        /* Grid overlap communication */
        if (npme[i] > 1)
        {
<<<<<<< HEAD
            nk = (i==0 ? ir->nkx : ir->nky);
            overlap = (nk % npme[i] == 0 ? ir->pme_order-1 : ir->pme_order);
            temp = npme[i];
            temp *= overlap;
            temp *= ir->nkx;
            temp *= ir->nky;
            temp *= ir->nkz;
            temp /= nk;
=======
            nk        = (i == 0 ? ir->nkx : ir->nky);
            overlap   = (nk % npme[i] == 0 ? ir->pme_order-1 : ir->pme_order);
            temp      = npme[i];
            temp     *= overlap;
            temp     *= ir->nkx;
            temp     *= ir->nky;
            temp     *= ir->nkz;
            temp     /= nk;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            comm_pme += temp;
/* Old line comm_pme += npme[i]*overlap*ir->nkx*ir->nky*ir->nkz/nk; */
        }
    }

    /* PME FFT communication volume.
     * This only takes the communication into account and not imbalance
     * in the calculation. But the imbalance in communication and calculation
     * are similar and therefore these formulas also prefer load balance
     * in the FFT and pme_solve calculation.
     */
    comm_pme += comm_pme_cost_vol(npme[YY], ir->nky, ir->nkz, ir->nkx);
    comm_pme += comm_pme_cost_vol(npme[XX], ir->nkx, ir->nky, ir->nkz);
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Add cost of pbc_dx for bondeds */
    cost_pbcdx = 0;
    if ((nc[XX] == 1 || nc[YY] == 1) || (nc[ZZ] == 1 && ir->ePBC != epbcXY))
    {
        if ((ddbox->tric_dir[XX] && nc[XX] == 1) ||
            (ddbox->tric_dir[YY] && nc[YY] == 1))
        {
            cost_pbcdx = pbcdxr*pbcdx_tric_fac;
        }
        else
        {
            cost_pbcdx = pbcdxr*pbcdx_rect_fac;
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (debug)
    {
        fprintf(debug,
                "nc %2d %2d %2d %2d %2d vol pp %6.4f pbcdx %6.4f pme %9.3e tot %9.3e\n",
<<<<<<< HEAD
                nc[XX],nc[YY],nc[ZZ],npme[XX],npme[YY],
                comm_vol,cost_pbcdx,comm_pme,
                3*natoms*(comm_vol + cost_pbcdx) + comm_pme);
    }
    
=======
                nc[XX], nc[YY], nc[ZZ], npme[XX], npme[YY],
                comm_vol, cost_pbcdx, comm_pme,
                3*natoms*(comm_vol + cost_pbcdx) + comm_pme);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return 3*natoms*(comm_vol + cost_pbcdx) + comm_pme;
}

static void assign_factors(gmx_domdec_t *dd,
<<<<<<< HEAD
                           real limit,real cutoff,
                           matrix box,gmx_ddbox_t *ddbox,
                           int natoms,t_inputrec *ir,
                           float pbcdxr,int npme,
                           int ndiv,int *div,int *mdiv,ivec ir_try,ivec opt)
{
    int x,y,z,i;
    float ce;
    
    if (ndiv == 0)
    {
        ce = comm_cost_est(dd,limit,cutoff,box,ddbox,
                           natoms,ir,pbcdxr,npme,ir_try);
        if (ce >= 0 && (opt[XX] == 0 ||
                        ce < comm_cost_est(dd,limit,cutoff,box,ddbox,
                                           natoms,ir,pbcdxr,
                                           npme,opt)))
        {
            copy_ivec(ir_try,opt);
        }
        
        return;
    }
    
    for(x=mdiv[0]; x>=0; x--)
    {
        for(i=0; i<x; i++)
        {
            ir_try[XX] *= div[0];
        }
        for(y=mdiv[0]-x; y>=0; y--)
        {
            for(i=0; i<y; i++)
            {
                ir_try[YY] *= div[0];
            }
            for(i=0; i<mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] *= div[0];
            }
            
            /* recurse */
            assign_factors(dd,limit,cutoff,box,ddbox,natoms,ir,pbcdxr,npme,
                           ndiv-1,div+1,mdiv+1,ir_try,opt);
            
            for(i=0; i<mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] /= div[0];
            }
            for(i=0; i<y; i++)
=======
                           real limit, real cutoff,
                           matrix box, gmx_ddbox_t *ddbox,
                           int natoms, t_inputrec *ir,
                           float pbcdxr, int npme,
                           int ndiv, int *div, int *mdiv, ivec ir_try, ivec opt)
{
    int   x, y, z, i;
    float ce;

    if (ndiv == 0)
    {
        ce = comm_cost_est(dd, limit, cutoff, box, ddbox,
                           natoms, ir, pbcdxr, npme, ir_try);
        if (ce >= 0 && (opt[XX] == 0 ||
                        ce < comm_cost_est(dd, limit, cutoff, box, ddbox,
                                           natoms, ir, pbcdxr,
                                           npme, opt)))
        {
            copy_ivec(ir_try, opt);
        }

        return;
    }

    for (x = mdiv[0]; x >= 0; x--)
    {
        for (i = 0; i < x; i++)
        {
            ir_try[XX] *= div[0];
        }
        for (y = mdiv[0]-x; y >= 0; y--)
        {
            for (i = 0; i < y; i++)
            {
                ir_try[YY] *= div[0];
            }
            for (i = 0; i < mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] *= div[0];
            }

            /* recurse */
            assign_factors(dd, limit, cutoff, box, ddbox, natoms, ir, pbcdxr, npme,
                           ndiv-1, div+1, mdiv+1, ir_try, opt);

            for (i = 0; i < mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] /= div[0];
            }
            for (i = 0; i < y; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                ir_try[YY] /= div[0];
            }
        }
<<<<<<< HEAD
        for(i=0; i<x; i++)
=======
        for (i = 0; i < x; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            ir_try[XX] /= div[0];
        }
    }
}

static real optimize_ncells(FILE *fplog,
<<<<<<< HEAD
                            int nnodes_tot,int npme_only,
                            gmx_bool bDynLoadBal,real dlb_scale,
                            gmx_mtop_t *mtop,matrix box,gmx_ddbox_t *ddbox,
                            t_inputrec *ir,
                            gmx_domdec_t *dd,
                            real cellsize_limit,real cutoff,
                            gmx_bool bInterCGBondeds,gmx_bool bInterCGMultiBody,
                            ivec nc)
{
    int npp,npme,ndiv,*div,*mdiv,d,nmax;
    gmx_bool bExcl_pbcdx;
    float pbcdxr;
    real limit;
    ivec itry;
    
    limit  = cellsize_limit;
    
=======
                            int nnodes_tot, int npme_only,
                            gmx_bool bDynLoadBal, real dlb_scale,
                            gmx_mtop_t *mtop, matrix box, gmx_ddbox_t *ddbox,
                            t_inputrec *ir,
                            gmx_domdec_t *dd,
                            real cellsize_limit, real cutoff,
                            gmx_bool bInterCGBondeds, gmx_bool bInterCGMultiBody,
                            ivec nc)
{
    int      npp, npme, ndiv, *div, *mdiv, d, nmax;
    gmx_bool bExcl_pbcdx;
    float    pbcdxr;
    real     limit;
    ivec     itry;

    limit  = cellsize_limit;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    dd->nc[XX] = 1;
    dd->nc[YY] = 1;
    dd->nc[ZZ] = 1;

    npp = nnodes_tot - npme_only;
    if (EEL_PME(ir->coulombtype))
    {
        npme = (npme_only > 0 ? npme_only : npp);
    }
    else
    {
        npme = 0;
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (bInterCGBondeds)
    {
        /* For Ewald exclusions pbc_dx is not called */
        bExcl_pbcdx =
            (IR_EXCL_FORCES(*ir) && !EEL_FULL(ir->coulombtype));
<<<<<<< HEAD
        pbcdxr = (double)n_bonded_dx(mtop,bExcl_pbcdx)/(double)mtop->natoms;
=======
        pbcdxr = (double)n_bonded_dx(mtop, bExcl_pbcdx)/(double)mtop->natoms;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        /* Every molecule is a single charge group: no pbc required */
        pbcdxr = 0;
    }
    /* Add a margin for DLB and/or pressure scaling */
    if (bDynLoadBal)
    {
        if (dlb_scale >= 1.0)
        {
<<<<<<< HEAD
            gmx_fatal(FARGS,"The value for option -dds should be smaller than 1");
        }
        if (fplog)
        {
            fprintf(fplog,"Scaling the initial minimum size with 1/%g (option -dds) = %g\n",dlb_scale,1/dlb_scale);
=======
            gmx_fatal(FARGS, "The value for option -dds should be smaller than 1");
        }
        if (fplog)
        {
            fprintf(fplog, "Scaling the initial minimum size with 1/%g (option -dds) = %g\n", dlb_scale, 1/dlb_scale);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        limit /= dlb_scale;
    }
    else if (ir->epc != epcNO)
    {
        if (fplog)
        {
<<<<<<< HEAD
            fprintf(fplog,"To account for pressure scaling, scaling the initial minimum size with %g\n",DD_GRID_MARGIN_PRES_SCALE);
            limit *= DD_GRID_MARGIN_PRES_SCALE;
        }
    }
    
    if (fplog)
    {
        fprintf(fplog,"Optimizing the DD grid for %d cells with a minimum initial size of %.3f nm\n",npp,limit);

        if (inhomogeneous_z(ir))
        {
            fprintf(fplog,"Ewald_geometry=%s: assuming inhomogeneous particle distribution in z, will not decompose in z.\n",eewg_names[ir->ewald_geometry]);
=======
            fprintf(fplog, "To account for pressure scaling, scaling the initial minimum size with %g\n", DD_GRID_MARGIN_PRES_SCALE);
            limit *= DD_GRID_MARGIN_PRES_SCALE;
        }
    }

    if (fplog)
    {
        fprintf(fplog, "Optimizing the DD grid for %d cells with a minimum initial size of %.3f nm\n", npp, limit);

        if (inhomogeneous_z(ir))
        {
            fprintf(fplog, "Ewald_geometry=%s: assuming inhomogeneous particle distribution in z, will not decompose in z.\n", eewg_names[ir->ewald_geometry]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }

        if (limit > 0)
        {
<<<<<<< HEAD
            fprintf(fplog,"The maximum allowed number of cells is:");
            for(d=0; d<DIM; d++)
=======
            fprintf(fplog, "The maximum allowed number of cells is:");
            for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                nmax = (int)(ddbox->box_size[d]*ddbox->skew_fac[d]/limit);
                if (d >= ddbox->npbcdim && nmax < 2)
                {
                    nmax = 2;
                }
                if (d == ZZ && inhomogeneous_z(ir))
                {
                    nmax = 1;
                }
<<<<<<< HEAD
                fprintf(fplog," %c %d",'X' + d,nmax);
            }
            fprintf(fplog,"\n");
        }
    }
    
    if (debug)
    {
        fprintf(debug,"Average nr of pbc_dx calls per atom %.2f\n",pbcdxr);
    }
    
    /* Decompose npp in factors */
    ndiv = factorize(npp,&div,&mdiv);
    
=======
                fprintf(fplog, " %c %d", 'X' + d, nmax);
            }
            fprintf(fplog, "\n");
        }
    }

    if (debug)
    {
        fprintf(debug, "Average nr of pbc_dx calls per atom %.2f\n", pbcdxr);
    }

    /* Decompose npp in factors */
    ndiv = factorize(npp, &div, &mdiv);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    itry[XX] = 1;
    itry[YY] = 1;
    itry[ZZ] = 1;
    clear_ivec(nc);
<<<<<<< HEAD
    assign_factors(dd,limit,cutoff,box,ddbox,mtop->natoms,ir,pbcdxr,
                   npme,ndiv,div,mdiv,itry,nc);
    
    sfree(div);
    sfree(mdiv);
    
=======
    assign_factors(dd, limit, cutoff, box, ddbox, mtop->natoms, ir, pbcdxr,
                   npme, ndiv, div, mdiv, itry, nc);

    sfree(div);
    sfree(mdiv);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return limit;
}

real dd_choose_grid(FILE *fplog,
<<<<<<< HEAD
                    t_commrec *cr,gmx_domdec_t *dd,t_inputrec *ir,
                    gmx_mtop_t *mtop,matrix box,gmx_ddbox_t *ddbox,
                    gmx_bool bDynLoadBal,real dlb_scale,
                    real cellsize_limit,real cutoff_dd,
                    gmx_bool bInterCGBondeds,gmx_bool bInterCGMultiBody)
{
    int  nnodes_div,ldiv;
    real limit;
    
=======
                    t_commrec *cr, gmx_domdec_t *dd, t_inputrec *ir,
                    gmx_mtop_t *mtop, matrix box, gmx_ddbox_t *ddbox,
                    gmx_bool bDynLoadBal, real dlb_scale,
                    real cellsize_limit, real cutoff_dd,
                    gmx_bool bInterCGBondeds, gmx_bool bInterCGMultiBody)
{
    gmx_large_int_t nnodes_div, ldiv;
    real            limit;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (MASTER(cr))
    {
        nnodes_div = cr->nnodes;
        if (EEL_PME(ir->coulombtype))
        {
            if (cr->npmenodes > 0)
            {
                if (cr->nnodes <= 2)
                {
                    gmx_fatal(FARGS,
                              "Can not have separate PME nodes with 2 or less nodes");
                }
                if (cr->npmenodes >= cr->nnodes)
                {
                    gmx_fatal(FARGS,
                              "Can not have %d separate PME nodes with just %d total nodes",
                              cr->npmenodes, cr->nnodes);
                }

                /* If the user purposely selected the number of PME nodes,
                 * only check for large primes in the PP node count.
                 */
                nnodes_div -= cr->npmenodes;
            }
        }
        else
        {
            cr->npmenodes = 0;
        }

        if (cr->nnodes > 12)
        {
            ldiv = largest_divisor(nnodes_div);
            /* Check if the largest divisor is more than nnodes^2/3 */
            if (ldiv*ldiv*ldiv > nnodes_div*nnodes_div)
            {
<<<<<<< HEAD
                gmx_fatal(FARGS,"The number of nodes you selected (%d) contains a large prime factor %d. In most cases this will lead to bad performance. Choose a number with smaller prime factors or set the decomposition (option -dd) manually.",
                          nnodes_div,ldiv);
=======
                gmx_fatal(FARGS, "The number of nodes you selected (%d) contains a large prime factor %d. In most cases this will lead to bad performance. Choose a number with smaller prime factors or set the decomposition (option -dd) manually.",
                          nnodes_div, ldiv);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }

        if (EEL_PME(ir->coulombtype))
        {
            if (cr->npmenodes < 0)
            {
<<<<<<< HEAD
                if (cr->nnodes <= 10)
                {
                    cr->npmenodes = 0;
                }
                else
                {
                    cr->npmenodes = guess_npme(fplog,mtop,ir,box,cr->nnodes);
                }
            }
            if (fplog)
            {
                fprintf(fplog,"Using %d separate PME nodes\n",cr->npmenodes);
            }
        }
        
        limit = optimize_ncells(fplog,cr->nnodes,cr->npmenodes,
                                bDynLoadBal,dlb_scale,
                                mtop,box,ddbox,ir,dd,
                                cellsize_limit,cutoff_dd,
                                bInterCGBondeds,bInterCGMultiBody,
=======
                /* Use PME nodes when the number of nodes is more than 16 */
                if (cr->nnodes <= 18)
                {
                    cr->npmenodes = 0;
                    if (fplog)
                    {
                        fprintf(fplog, "Using %d separate PME nodes, as there are too few total\n nodes for efficient splitting\n", cr->npmenodes);
                    }
                }
                else
                {
                    cr->npmenodes = guess_npme(fplog, mtop, ir, box, cr->nnodes);
                    if (fplog)
                    {
                        fprintf(fplog, "Using %d separate PME nodes, as guessed by mdrun\n", cr->npmenodes);
                    }
                }
            }
            else
            {
                if (fplog)
                {
                    fprintf(fplog, "Using %d separate PME nodes, per user request\n", cr->npmenodes);
                }
            }
        }

        limit = optimize_ncells(fplog, cr->nnodes, cr->npmenodes,
                                bDynLoadBal, dlb_scale,
                                mtop, box, ddbox, ir, dd,
                                cellsize_limit, cutoff_dd,
                                bInterCGBondeds, bInterCGMultiBody,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                dd->nc);
    }
    else
    {
        limit = 0;
    }
    /* Communicate the information set by the master to all nodes */
<<<<<<< HEAD
    gmx_bcast(sizeof(dd->nc),dd->nc,cr);
    if (EEL_PME(ir->coulombtype))
    {
        gmx_bcast(sizeof(ir->nkx),&ir->nkx,cr);
        gmx_bcast(sizeof(ir->nky),&ir->nky,cr);
        gmx_bcast(sizeof(cr->npmenodes),&cr->npmenodes,cr);
=======
    gmx_bcast(sizeof(dd->nc), dd->nc, cr);
    if (EEL_PME(ir->coulombtype))
    {
        gmx_bcast(sizeof(ir->nkx), &ir->nkx, cr);
        gmx_bcast(sizeof(ir->nky), &ir->nky, cr);
        gmx_bcast(sizeof(cr->npmenodes), &cr->npmenodes, cr);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        cr->npmenodes = 0;
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return limit;
}
