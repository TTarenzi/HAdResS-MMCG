<<<<<<< HEAD
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
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
 * GROwing Monsters And Cloning Shrimps
=======
/*
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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <math.h>
#include "macros.h"
#include "main.h"
#include "smalloc.h"
#include "futil.h"
#include "tgroup.h"
#include "vec.h"
#include "network.h"
#include "smalloc.h"
#include "update.h"
#include "rbin.h"
#include "mtop_util.h"
<<<<<<< HEAD

static void init_grptcstat(int ngtc,t_grp_tcstat tcstat[])
{ 
    int i,j;
    
    for(i=0; (i<ngtc); i++) {
=======
#include "gmx_omp_nthreads.h"

static void init_grptcstat(int ngtc, t_grp_tcstat tcstat[])
{
    int i, j;

    for (i = 0; (i < ngtc); i++)
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        tcstat[i].T = 0;
        clear_mat(tcstat[i].ekinh);
        clear_mat(tcstat[i].ekinh_old);
        clear_mat(tcstat[i].ekinf);
    }
}

static void init_grpstat(FILE *log,
<<<<<<< HEAD
			 gmx_mtop_t *mtop,int ngacc,t_grp_acc gstat[])
{
    gmx_groups_t *groups;
    gmx_mtop_atomloop_all_t aloop;
    int    i,grp;
    t_atom *atom;
    
    if (ngacc > 0) {
        groups = &mtop->groups;
        aloop = gmx_mtop_atomloop_all_init(mtop);
        while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) 
        {
            grp = ggrpnr(groups,egcACC,i);
=======
                         gmx_mtop_t *mtop, int ngacc, t_grp_acc gstat[])
{
    gmx_groups_t           *groups;
    gmx_mtop_atomloop_all_t aloop;
    int                     i, grp;
    t_atom                 *atom;

    if (ngacc > 0)
    {
        groups = &mtop->groups;
        aloop  = gmx_mtop_atomloop_all_init(mtop);
        while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
        {
            grp = ggrpnr(groups, egcACC, i);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if ((grp < 0) && (grp >= ngacc))
            {
                gmx_incons("Input for acceleration groups wrong");
            }
            gstat[grp].nat++;
            /* This will not work for integrator BD */
            gstat[grp].mA += atom->m;
            gstat[grp].mB += atom->mB;
        }
    }
}

<<<<<<< HEAD
void init_ekindata(FILE *log,gmx_mtop_t *mtop,t_grpopts *opts,
                   gmx_ekindata_t *ekind)
{
  int i;
#ifdef DEBUG
  fprintf(log,"ngtc: %d, ngacc: %d, ngener: %d\n",opts->ngtc,opts->ngacc,
	  opts->ngener);
#endif

  /* bNEMD tells if we should remove remove the COM velocity
   * from the velocities during velocity scaling in T-coupling.
   * Turn this on when we have multiple acceleration groups
   * or one accelerated group.
   */
  ekind->bNEMD = (opts->ngacc > 1 || norm(opts->acc[0]) > 0);

  ekind->ngtc = opts->ngtc;
  snew(ekind->tcstat,opts->ngtc);
  init_grptcstat(opts->ngtc,ekind->tcstat);
  /* Set Berendsen tcoupl lambda's to 1, 
   * so runs without Berendsen coupling are not affected.
   */
  for(i=0; i<opts->ngtc; i++) 
  {
      ekind->tcstat[i].lambda = 1.0;
      ekind->tcstat[i].vscale_nhc = 1.0;
      ekind->tcstat[i].ekinscaleh_nhc = 1.0;
      ekind->tcstat[i].ekinscalef_nhc = 1.0;
  }
  
  ekind->ngacc = opts->ngacc;
  snew(ekind->grpstat,opts->ngacc);
  init_grpstat(log,mtop,opts->ngacc,ekind->grpstat);
}

void accumulate_u(t_commrec *cr,t_grpopts *opts,gmx_ekindata_t *ekind)
{
    /* This routine will only be called when it's necessary */
    t_bin *rb;
    int   g;
    
    rb = mk_bin();
    
    for(g=0; (g<opts->ngacc); g++) 
    {
        add_binr(rb,DIM,ekind->grpstat[g].u);
    }
    sum_bin(rb,cr);
    
    for(g=0; (g<opts->ngacc); g++) 
    {
        extract_binr(rb,DIM*g,DIM,ekind->grpstat[g].u);
    }
    destroy_bin(rb);
}       
=======
void init_ekindata(FILE *log, gmx_mtop_t *mtop, t_grpopts *opts,
                   gmx_ekindata_t *ekind)
{
    int i;
    int nthread, thread;
#ifdef DEBUG
    fprintf(log, "ngtc: %d, ngacc: %d, ngener: %d\n", opts->ngtc, opts->ngacc,
            opts->ngener);
#endif

    /* bNEMD tells if we should remove remove the COM velocity
     * from the velocities during velocity scaling in T-coupling.
     * Turn this on when we have multiple acceleration groups
     * or one accelerated group.
     */
    ekind->bNEMD = (opts->ngacc > 1 || norm(opts->acc[0]) > 0);

    ekind->ngtc = opts->ngtc;
    snew(ekind->tcstat, opts->ngtc);
    init_grptcstat(opts->ngtc, ekind->tcstat);
    /* Set Berendsen tcoupl lambda's to 1,
     * so runs without Berendsen coupling are not affected.
     */
    for (i = 0; i < opts->ngtc; i++)
    {
        ekind->tcstat[i].lambda         = 1.0;
        ekind->tcstat[i].vscale_nhc     = 1.0;
        ekind->tcstat[i].ekinscaleh_nhc = 1.0;
        ekind->tcstat[i].ekinscalef_nhc = 1.0;
    }

    nthread = gmx_omp_nthreads_get(emntUpdate);

    snew(ekind->ekin_work_alloc, nthread);
    snew(ekind->ekin_work, nthread);
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        /* Allocate 2 elements extra on both sides,
         * so in single precision we have 2*3*3*4=72 bytes buffer
         * on both sides to avoid cache pollution.
         */
        snew(ekind->ekin_work_alloc[thread], ekind->ngtc+4);
        ekind->ekin_work[thread] = ekind->ekin_work_alloc[thread] + 2;
    }

    ekind->ngacc = opts->ngacc;
    snew(ekind->grpstat, opts->ngacc);
    init_grpstat(log, mtop, opts->ngacc, ekind->grpstat);
}

void accumulate_u(t_commrec *cr, t_grpopts *opts, gmx_ekindata_t *ekind)
{
    /* This routine will only be called when it's necessary */
    t_bin *rb;
    int    g;

    rb = mk_bin();

    for (g = 0; (g < opts->ngacc); g++)
    {
        add_binr(rb, DIM, ekind->grpstat[g].u);
    }
    sum_bin(rb, cr);

    for (g = 0; (g < opts->ngacc); g++)
    {
        extract_binr(rb, DIM*g, DIM, ekind->grpstat[g].u);
    }
    destroy_bin(rb);
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/* I don't think accumulate_ekin is used anymore? */

#if 0
<<<<<<< HEAD
static void accumulate_ekin(t_commrec *cr,t_grpopts *opts,
			    gmx_ekindata_t *ekind)
{
    int g;

    if(PAR(cr))
    {
        for(g=0; (g<opts->ngtc); g++) 
        {
            gmx_sum(DIM*DIM,ekind->tcstat[g].ekinf[0],cr);
        }
    }
}       
#endif 

void update_ekindata(int start,int homenr,gmx_ekindata_t *ekind,
                     t_grpopts *opts,rvec v[],t_mdatoms *md,real lambda)
{
  int  d,g,n;
  real mv;

  /* calculate mean velocities at whole timestep */ 
  for(g=0; (g<opts->ngtc); g++) {
    ekind->tcstat[g].T = 0;
  }

  if (ekind->bNEMD) {
    for (g=0; (g<opts->ngacc); g++)
      clear_rvec(ekind->grpstat[g].u);
    
    g = 0;
    for(n=start; (n<start+homenr); n++) {
      if (md->cACC)
	g = md->cACC[n];
      for(d=0; (d<DIM);d++) {
          mv = md->massT[n]*v[n][d];
          ekind->grpstat[g].u[d] += mv;
      }
    }

    for (g=0; (g < opts->ngacc); g++) {
      for(d=0; (d<DIM);d++) {
          ekind->grpstat[g].u[d] /=
              (1-lambda)*ekind->grpstat[g].mA + lambda*ekind->grpstat[g].mB;
      }
    }
  }
}

real sum_ekin(t_grpopts *opts,gmx_ekindata_t *ekind,real *dekindlambda,
              gmx_bool bEkinAveVel, gmx_bool bSaveEkinOld, gmx_bool bScaleEkin)
{
    int          i,j,m,ngtc;
    real         T,ek;
    t_grp_tcstat *tcstat;
    real         nrdf,nd,*ndf;
    
    ngtc = opts->ngtc;
    ndf  = opts->nrdf;
    
    T = 0; 
=======
static void accumulate_ekin(t_commrec *cr, t_grpopts *opts,
                            gmx_ekindata_t *ekind)
{
    int g;

    if (PAR(cr))
    {
        for (g = 0; (g < opts->ngtc); g++)
        {
            gmx_sum(DIM*DIM, ekind->tcstat[g].ekinf[0], cr);
        }
    }
}
#endif

void update_ekindata(int start, int homenr, gmx_ekindata_t *ekind,
                     t_grpopts *opts, rvec v[], t_mdatoms *md, real lambda)
{
    int  d, g, n;
    real mv;

    /* calculate mean velocities at whole timestep */
    for (g = 0; (g < opts->ngtc); g++)
    {
        ekind->tcstat[g].T = 0;
    }

    if (ekind->bNEMD)
    {
        for (g = 0; (g < opts->ngacc); g++)
        {
            clear_rvec(ekind->grpstat[g].u);
        }

        g = 0;
        for (n = start; (n < start+homenr); n++)
        {
            if (md->cACC)
            {
                g = md->cACC[n];
            }
            for (d = 0; (d < DIM); d++)
            {
                mv                      = md->massT[n]*v[n][d];
                ekind->grpstat[g].u[d] += mv;
            }
        }

        for (g = 0; (g < opts->ngacc); g++)
        {
            for (d = 0; (d < DIM); d++)
            {
                ekind->grpstat[g].u[d] /=
                    (1-lambda)*ekind->grpstat[g].mA + lambda*ekind->grpstat[g].mB;
            }
        }
    }
}

real sum_ekin(t_grpopts *opts, gmx_ekindata_t *ekind, real *dekindlambda,
              gmx_bool bEkinAveVel, gmx_bool bSaveEkinOld, gmx_bool bScaleEkin)
{
    int           i, j, m, ngtc;
    real          T, ek;
    t_grp_tcstat *tcstat;
    real          nrdf, nd, *ndf;

    ngtc = opts->ngtc;
    ndf  = opts->nrdf;

    T    = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    nrdf = 0;

    clear_mat(ekind->ekin);

<<<<<<< HEAD
    for(i=0; (i<ngtc); i++) 
    {

        nd = ndf[i];
=======
    for (i = 0; (i < ngtc); i++)
    {

        nd     = ndf[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        tcstat = &ekind->tcstat[i];
        /* Sometimes a group does not have degrees of freedom, e.g.
         * when it consists of shells and virtual sites, then we just
         * set the temperatue to 0 and also neglect the kinetic
         * energy, which should be  zero anyway.
         */
<<<<<<< HEAD
        
        if (nd > 0) {
=======

        if (nd > 0)
        {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (bEkinAveVel)
            {
                if (!bScaleEkin)
                {
                    /* in this case, kinetic energy is from the current velocities already */
<<<<<<< HEAD
                    msmul(tcstat->ekinf,tcstat->ekinscalef_nhc,tcstat->ekinf);
                }
            } 
            else 
                
            {
                /* Calculate the full step Ekin as the average of the half steps */
                for(j=0; (j<DIM); j++)
                {
                    for(m=0; (m<DIM); m++)
=======
                    msmul(tcstat->ekinf, tcstat->ekinscalef_nhc, tcstat->ekinf);
                }
            }
            else

            {
                /* Calculate the full step Ekin as the average of the half steps */
                for (j = 0; (j < DIM); j++)
                {
                    for (m = 0; (m < DIM); m++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        tcstat->ekinf[j][m] =
                            0.5*(tcstat->ekinh[j][m]*tcstat->ekinscaleh_nhc + tcstat->ekinh_old[j][m]);
                    }
                }
            }
<<<<<<< HEAD
            m_add(tcstat->ekinf,ekind->ekin,ekind->ekin);
            
            tcstat->Th = calc_temp(trace(tcstat->ekinh),nd);
            tcstat->T  = calc_temp(trace(tcstat->ekinf),nd);

            /* after the scaling factors have been multiplied in, we can remove them */
            if (bEkinAveVel) 
            {
                tcstat->ekinscalef_nhc = 1.0;
            } 
            else 
=======
            m_add(tcstat->ekinf, ekind->ekin, ekind->ekin);

            tcstat->Th = calc_temp(trace(tcstat->ekinh), nd);
            tcstat->T  = calc_temp(trace(tcstat->ekinf), nd);

            /* after the scaling factors have been multiplied in, we can remove them */
            if (bEkinAveVel)
            {
                tcstat->ekinscalef_nhc = 1.0;
            }
            else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                tcstat->ekinscaleh_nhc = 1.0;
            }
        }
<<<<<<< HEAD
        else 
=======
        else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            tcstat->T  = 0;
            tcstat->Th = 0;
        }
        T    += nd*tcstat->T;
        nrdf += nd;
    }
    if (nrdf > 0)
    {
<<<<<<< HEAD
        T/=nrdf;
    }
    if (dekindlambda) 
=======
        T /= nrdf;
    }
    if (dekindlambda)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        *dekindlambda = 0.5*(ekind->dekindl + ekind->dekindl_old);
    }
    return T;
}