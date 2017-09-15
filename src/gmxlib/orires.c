<<<<<<< HEAD
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-  
 *
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
 * GROningen Mixture of Alchemy and Childrens' Stories
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "nrjac.h"
#include "network.h"
#include "orires.h"
#include "do_fit.h"
#include "main.h"
#include "copyrite.h"
#include "pbc.h"
#include "mtop_util.h"

<<<<<<< HEAD
void init_orires(FILE *fplog,const gmx_mtop_t *mtop,
                 rvec xref[],
                 const t_inputrec *ir,
                 const gmx_multisim_t *ms,t_oriresdata *od,
                 t_state *state)
{
    int    i,j,d,ex,nmol,nr,*nr_ex;
    double mtot;
    rvec   com;
    gmx_mtop_ilistloop_t iloop;
    t_ilist *il;
    gmx_mtop_atomloop_all_t aloop;
    t_atom *atom;
=======
void init_orires(FILE *fplog, const gmx_mtop_t *mtop,
                 rvec xref[],
                 const t_inputrec *ir,
                 const gmx_multisim_t *ms, t_oriresdata *od,
                 t_state *state)
{
    int                     i, j, d, ex, nmol, nr, *nr_ex;
    double                  mtot;
    rvec                    com;
    gmx_mtop_ilistloop_t    iloop;
    t_ilist                *il;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    od->fc  = ir->orires_fc;
    od->nex = 0;
    od->S   = NULL;

<<<<<<< HEAD
    od->M=NULL;
    od->eig=NULL;
    od->v=NULL;

    od->nr = gmx_mtop_ftype_count(mtop,F_ORIRES);
=======
    od->M   = NULL;
    od->eig = NULL;
    od->v   = NULL;

    od->nr = gmx_mtop_ftype_count(mtop, F_ORIRES);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (od->nr == 0)
    {
        return;
    }
<<<<<<< HEAD
    
    nr_ex = NULL;
    
    iloop = gmx_mtop_ilistloop_init(mtop);
    while (gmx_mtop_ilistloop_next(iloop,&il,&nmol))
    {
        for(i=0; i<il[F_ORIRES].nr; i+=3)
=======

    nr_ex = NULL;

    iloop = gmx_mtop_ilistloop_init(mtop);
    while (gmx_mtop_ilistloop_next(iloop, &il, &nmol))
    {
        for (i = 0; i < il[F_ORIRES].nr; i += 3)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            ex = mtop->ffparams.iparams[il[F_ORIRES].iatoms[i]].orires.ex;
            if (ex >= od->nex)
            {
<<<<<<< HEAD
                srenew(nr_ex,ex+1);
                for(j=od->nex; j<ex+1; j++)
                {
                    nr_ex[j] = 0;
            }
=======
                srenew(nr_ex, ex+1);
                for (j = od->nex; j < ex+1; j++)
                {
                    nr_ex[j] = 0;
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                od->nex = ex+1;
            }
            nr_ex[ex]++;
        }
    }
<<<<<<< HEAD
    snew(od->S,od->nex);
    /* When not doing time averaging, the instaneous and time averaged data
     * are indentical and the pointers can point to the same memory.
     */
    snew(od->Dinsl,od->nr);
    if (ms)
    {
        snew(od->Dins,od->nr);
=======
    snew(od->S, od->nex);
    /* When not doing time averaging, the instaneous and time averaged data
     * are indentical and the pointers can point to the same memory.
     */
    snew(od->Dinsl, od->nr);
    if (ms)
    {
        snew(od->Dins, od->nr);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        od->Dins = od->Dinsl;
    }

    if (ir->orires_tau == 0)
    {
<<<<<<< HEAD
        od->Dtav = od->Dins;
        od->edt  = 0.0;
        od->edt_1= 1.0;
    }
    else
    {
        snew(od->Dtav,od->nr);
        od->edt  = exp(-ir->delta_t/ir->orires_tau);
        od->edt_1= 1.0 - od->edt;

        /* Extend the state with the orires history */
        state->flags |= (1<<estORIRE_INITF);
        state->hist.orire_initf = 1;
        state->flags |= (1<<estORIRE_DTAV);
        state->hist.norire_Dtav = od->nr*5;
        snew(state->hist.orire_Dtav,state->hist.norire_Dtav);
    }

    snew(od->oinsl,od->nr);
    if (ms)
    {
        snew(od->oins,od->nr);
=======
        od->Dtav  = od->Dins;
        od->edt   = 0.0;
        od->edt_1 = 1.0;
    }
    else
    {
        snew(od->Dtav, od->nr);
        od->edt   = exp(-ir->delta_t/ir->orires_tau);
        od->edt_1 = 1.0 - od->edt;

        /* Extend the state with the orires history */
        state->flags           |= (1<<estORIRE_INITF);
        state->hist.orire_initf = 1;
        state->flags           |= (1<<estORIRE_DTAV);
        state->hist.norire_Dtav = od->nr*5;
        snew(state->hist.orire_Dtav, state->hist.norire_Dtav);
    }

    snew(od->oinsl, od->nr);
    if (ms)
    {
        snew(od->oins, od->nr);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        od->oins = od->oinsl;
    }
<<<<<<< HEAD
    if (ir->orires_tau == 0) {
=======
    if (ir->orires_tau == 0)
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        od->otav = od->oins;
    }
    else
    {
<<<<<<< HEAD
        snew(od->otav,od->nr);
    }
    snew(od->tmp,od->nex);
    snew(od->TMP,od->nex);
    for(ex=0; ex<od->nex; ex++)
    {
        snew(od->TMP[ex],5);
        for(i=0; i<5; i++)
        {
            snew(od->TMP[ex][i],5);
        }
    }
    
    od->nref = 0;
    for(i=0; i<mtop->natoms; i++)
    {
        if (ggrpnr(&mtop->groups,egcORFIT,i) == 0)
=======
        snew(od->otav, od->nr);
    }
    snew(od->tmp, od->nex);
    snew(od->TMP, od->nex);
    for (ex = 0; ex < od->nex; ex++)
    {
        snew(od->TMP[ex], 5);
        for (i = 0; i < 5; i++)
        {
            snew(od->TMP[ex][i], 5);
        }
    }

    od->nref = 0;
    for (i = 0; i < mtop->natoms; i++)
    {
        if (ggrpnr(&mtop->groups, egcORFIT, i) == 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            od->nref++;
        }
    }
<<<<<<< HEAD
    snew(od->mref,od->nref);
    snew(od->xref,od->nref);
    snew(od->xtmp,od->nref);
    
    snew(od->eig,od->nex*12);
    
=======
    snew(od->mref, od->nref);
    snew(od->xref, od->nref);
    snew(od->xtmp, od->nref);

    snew(od->eig, od->nex*12);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Determine the reference structure on the master node.
     * Copy it to the other nodes after checking multi compatibility,
     * so we are sure the subsystems match before copying.
     */
    clear_rvec(com);
<<<<<<< HEAD
    mtot = 0.0;
    j = 0;
    aloop = gmx_mtop_atomloop_all_init(mtop);
    while(gmx_mtop_atomloop_all_next(aloop,&i,&atom))
=======
    mtot  = 0.0;
    j     = 0;
    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (mtop->groups.grpnr[egcORFIT] == NULL ||
            mtop->groups.grpnr[egcORFIT][i] == 0)
        {
            /* Not correct for free-energy with changing masses */
            od->mref[j] = atom->m;
<<<<<<< HEAD
            if (ms==NULL || MASTERSIM(ms))
            {
                copy_rvec(xref[i],od->xref[j]);
                for(d=0; d<DIM; d++)
=======
            if (ms == NULL || MASTERSIM(ms))
            {
                copy_rvec(xref[i], od->xref[j]);
                for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    com[d] += od->mref[j]*xref[i][d];
                }
            }
            mtot += od->mref[j];
            j++;
        }
    }
<<<<<<< HEAD
    svmul(1.0/mtot,com,com);
    if (ms==NULL || MASTERSIM(ms))
    {
        for(j=0; j<od->nref; j++)
        {
            rvec_dec(od->xref[j],com);
        }
    }
    
    fprintf(fplog,"Found %d orientation experiments\n",od->nex);
    for(i=0; i<od->nex; i++)
    {
        fprintf(fplog,"  experiment %d has %d restraints\n",i+1,nr_ex[i]);
    }
    
    sfree(nr_ex);
    
    fprintf(fplog,"  the fit group consists of %d atoms and has total mass %g\n",
            od->nref,mtot);
    
    if (ms)
    {
        fprintf(fplog,"  the orientation restraints are ensemble averaged over %d systems\n",ms->nsim);
        
        check_multi_int(fplog,ms,od->nr,
                        "the number of orientation restraints");
        check_multi_int(fplog,ms,od->nref,
                        "the number of fit atoms for orientation restraining");
        check_multi_int(fplog,ms,ir->nsteps,"nsteps");
        /* Copy the reference coordinates from the master to the other nodes */
        gmx_sum_sim(DIM*od->nref,od->xref[0],ms);
    }
    
    please_cite(fplog,"Hess2003");
=======
    svmul(1.0/mtot, com, com);
    if (ms == NULL || MASTERSIM(ms))
    {
        for (j = 0; j < od->nref; j++)
        {
            rvec_dec(od->xref[j], com);
        }
    }

    fprintf(fplog, "Found %d orientation experiments\n", od->nex);
    for (i = 0; i < od->nex; i++)
    {
        fprintf(fplog, "  experiment %d has %d restraints\n", i+1, nr_ex[i]);
    }

    sfree(nr_ex);

    fprintf(fplog, "  the fit group consists of %d atoms and has total mass %g\n",
            od->nref, mtot);

    if (ms)
    {
        fprintf(fplog, "  the orientation restraints are ensemble averaged over %d systems\n", ms->nsim);

        check_multi_int(fplog, ms, od->nr,
                        "the number of orientation restraints",
                        FALSE);
        check_multi_int(fplog, ms, od->nref,
                        "the number of fit atoms for orientation restraining",
                        FALSE);
        check_multi_int(fplog, ms, ir->nsteps, "nsteps", FALSE);
        /* Copy the reference coordinates from the master to the other nodes */
        gmx_sum_sim(DIM*od->nref, od->xref[0], ms);
    }

    please_cite(fplog, "Hess2003");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

void diagonalize_orires_tensors(t_oriresdata *od)
{
<<<<<<< HEAD
    int           ex,i,j,nrot,ord[DIM],t;
    matrix        S,TMP;
    
    if (od->M == NULL)
    {
        snew(od->M,DIM);
        for(i=0; i<DIM; i++)
        {
            snew(od->M[i],DIM);
        }
        snew(od->eig_diag,DIM);
        snew(od->v,DIM);
        for(i=0; i<DIM; i++)
        {
            snew(od->v[i],DIM);
        }
    }

    for(ex=0; ex<od->nex; ex++)
    {
        /* Rotate the S tensor back to the reference frame */
        mmul(od->R,od->S[ex],TMP);
        mtmul(TMP,od->R,S);
        for(i=0; i<DIM; i++)
        {
            for(j=0; j<DIM; j++)
=======
    int           ex, i, j, nrot, ord[DIM], t;
    matrix        S, TMP;

    if (od->M == NULL)
    {
        snew(od->M, DIM);
        for (i = 0; i < DIM; i++)
        {
            snew(od->M[i], DIM);
        }
        snew(od->eig_diag, DIM);
        snew(od->v, DIM);
        for (i = 0; i < DIM; i++)
        {
            snew(od->v[i], DIM);
        }
    }

    for (ex = 0; ex < od->nex; ex++)
    {
        /* Rotate the S tensor back to the reference frame */
        mmul(od->R, od->S[ex], TMP);
        mtmul(TMP, od->R, S);
        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                od->M[i][j] = S[i][j];
            }
        }
<<<<<<< HEAD
        
        jacobi(od->M,DIM,od->eig_diag,od->v,&nrot);
        
        for(i=0; i<DIM; i++)
        {
            ord[i] = i;
        }
        for(i=0; i<DIM; i++)
        {
            for(j=i+1; j<DIM; j++)
            {
                if (sqr(od->eig_diag[ord[j]]) > sqr(od->eig_diag[ord[i]]))
                {
                    t = ord[i];
=======

        jacobi(od->M, DIM, od->eig_diag, od->v, &nrot);

        for (i = 0; i < DIM; i++)
        {
            ord[i] = i;
        }
        for (i = 0; i < DIM; i++)
        {
            for (j = i+1; j < DIM; j++)
            {
                if (sqr(od->eig_diag[ord[j]]) > sqr(od->eig_diag[ord[i]]))
                {
                    t      = ord[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    ord[i] = ord[j];
                    ord[j] = t;
                }
            }
        }
<<<<<<< HEAD
            
        for(i=0; i<DIM; i++)
        {
            od->eig[ex*12 + i] = od->eig_diag[ord[i]];
        }
        for(i=0; i<DIM; i++)
        {
            for(j=0; j<DIM; j++)
=======

        for (i = 0; i < DIM; i++)
        {
            od->eig[ex*12 + i] = od->eig_diag[ord[i]];
        }
        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                od->eig[ex*12 + 3 + 3*i + j] = od->v[j][ord[i]];
            }
        }
    }
}

<<<<<<< HEAD
void print_orires_log(FILE *log,t_oriresdata *od)
{
    int  ex,i;
    real *eig;      
    
    diagonalize_orires_tensors(od);
    
    for(ex=0; ex<od->nex; ex++)
    {
        eig = od->eig + ex*12;
        fprintf(log,"  Orientation experiment %d:\n",ex+1);
        fprintf(log,"    order parameter: %g\n",eig[0]);
        for(i=0; i<DIM; i++)
        {
            fprintf(log,"    eig: %6.3f   %6.3f %6.3f %6.3f\n",
=======
void print_orires_log(FILE *log, t_oriresdata *od)
{
    int   ex, i;
    real *eig;

    diagonalize_orires_tensors(od);

    for (ex = 0; ex < od->nex; ex++)
    {
        eig = od->eig + ex*12;
        fprintf(log, "  Orientation experiment %d:\n", ex+1);
        fprintf(log, "    order parameter: %g\n", eig[0]);
        for (i = 0; i < DIM; i++)
        {
            fprintf(log, "    eig: %6.3f   %6.3f %6.3f %6.3f\n",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    (eig[0] != 0) ? eig[i]/eig[0] : eig[i],
                    eig[DIM+i*DIM+XX],
                    eig[DIM+i*DIM+YY],
                    eig[DIM+i*DIM+ZZ]);
        }
<<<<<<< HEAD
        fprintf(log,"\n");
=======
        fprintf(log, "\n");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}

real calc_orires_dev(const gmx_multisim_t *ms,
<<<<<<< HEAD
                     int nfa,const t_iatom forceatoms[],const t_iparams ip[],
                     const t_mdatoms *md,const rvec x[],const t_pbc *pbc,
                     t_fcdata *fcd,history_t *hist)
{
    int          fa,d,i,j,type,ex,nref;
    real         edt,edt_1,invn,pfac,r2,invr,corrfac,weight,wsv2,sw,dev;
    tensor       *S,R,TMP;
    rvec5        *Dinsl,*Dins,*Dtav,*rhs;
    real         *mref,***T;
    double       mtot;
    rvec         *xref,*xtmp,com,r_unrot,r;
    t_oriresdata *od;
    gmx_bool         bTAV;
    const real   two_thr=2.0/3.0;
    
=======
                     int nfa, const t_iatom forceatoms[], const t_iparams ip[],
                     const t_mdatoms *md, const rvec x[], const t_pbc *pbc,
                     t_fcdata *fcd, history_t *hist)
{
    int              fa, d, i, j, type, ex, nref;
    real             edt, edt_1, invn, pfac, r2, invr, corrfac, weight, wsv2, sw, dev;
    tensor          *S, R, TMP;
    rvec5           *Dinsl, *Dins, *Dtav, *rhs;
    real            *mref, ***T;
    double           mtot;
    rvec            *xref, *xtmp, com, r_unrot, r;
    t_oriresdata    *od;
    gmx_bool         bTAV;
    const real       two_thr = 2.0/3.0;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    od = &(fcd->orires);

    if (od->nr == 0)
    {
        /* This means that this is not the master node */
<<<<<<< HEAD
        gmx_fatal(FARGS,"Orientation restraints are only supported on the master node, use less processors");
    }
    
    bTAV = (od->edt != 0);
    edt  = od->edt;
    edt_1= od->edt_1;
    S    = od->S;
    Dinsl= od->Dinsl;
    Dins = od->Dins;
    Dtav = od->Dtav;
    T    = od->TMP;
    rhs  = od->tmp;
    nref = od->nref;
    mref = od->mref;
    xref = od->xref;
    xtmp = od->xtmp;
    
    if (bTAV)
    {
        od->exp_min_t_tau = hist->orire_initf*edt;
        
=======
        gmx_fatal(FARGS, "Orientation restraints are only supported on the master node, use less processors");
    }

    bTAV  = (od->edt != 0);
    edt   = od->edt;
    edt_1 = od->edt_1;
    S     = od->S;
    Dinsl = od->Dinsl;
    Dins  = od->Dins;
    Dtav  = od->Dtav;
    T     = od->TMP;
    rhs   = od->tmp;
    nref  = od->nref;
    mref  = od->mref;
    xref  = od->xref;
    xtmp  = od->xtmp;

    if (bTAV)
    {
        od->exp_min_t_tau = hist->orire_initf*edt;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Correction factor to correct for the lack of history
         * at short times.
         */
        corrfac = 1.0/(1.0 - od->exp_min_t_tau);
    }
    else
    {
        corrfac = 1.0;
    }

    if (ms)
    {
        invn = 1.0/ms->nsim;
    }
    else
    {
        invn = 1.0;
    }
<<<<<<< HEAD
    
    clear_rvec(com);
    mtot = 0;
    j=0;
    for(i=0; i<md->nr; i++)
    {
        if (md->cORF[i] == 0)
        {
            copy_rvec(x[i],xtmp[j]);
            mref[j] = md->massT[i];
            for(d=0; d<DIM; d++)
=======

    clear_rvec(com);
    mtot = 0;
    j    = 0;
    for (i = 0; i < md->nr; i++)
    {
        if (md->cORF[i] == 0)
        {
            copy_rvec(x[i], xtmp[j]);
            mref[j] = md->massT[i];
            for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                com[d] += mref[j]*xref[j][d];
            }
            mtot += mref[j];
            j++;
        }
    }
<<<<<<< HEAD
    svmul(1.0/mtot,com,com);
    for(j=0; j<nref; j++)
    {
        rvec_dec(xtmp[j],com);
    }
    /* Calculate the rotation matrix to rotate x to the reference orientation */
    calc_fit_R(DIM,nref,mref,xref,xtmp,R);
    copy_mat(R,od->R);
    
    d = 0;
    for(fa=0; fa<nfa; fa+=3)
=======
    svmul(1.0/mtot, com, com);
    for (j = 0; j < nref; j++)
    {
        rvec_dec(xtmp[j], com);
    }
    /* Calculate the rotation matrix to rotate x to the reference orientation */
    calc_fit_R(DIM, nref, mref, xref, xtmp, R);
    copy_mat(R, od->R);

    d = 0;
    for (fa = 0; fa < nfa; fa += 3)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        type = forceatoms[fa];
        if (pbc)
        {
<<<<<<< HEAD
            pbc_dx_aiuc(pbc,x[forceatoms[fa+1]],x[forceatoms[fa+2]],r_unrot);
        }
        else
        {
            rvec_sub(x[forceatoms[fa+1]],x[forceatoms[fa+2]],r_unrot);
        }
        mvmul(R,r_unrot,r);
=======
            pbc_dx_aiuc(pbc, x[forceatoms[fa+1]], x[forceatoms[fa+2]], r_unrot);
        }
        else
        {
            rvec_sub(x[forceatoms[fa+1]], x[forceatoms[fa+2]], r_unrot);
        }
        mvmul(R, r_unrot, r);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        r2   = norm2(r);
        invr = gmx_invsqrt(r2);
        /* Calculate the prefactor for the D tensor, this includes the factor 3! */
        pfac = ip[type].orires.c*invr*invr*3;
<<<<<<< HEAD
        for(i=0; i<ip[type].orires.power; i++)
=======
        for (i = 0; i < ip[type].orires.power; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            pfac *= invr;
        }
        Dinsl[d][0] = pfac*(2*r[0]*r[0] + r[1]*r[1] - r2);
        Dinsl[d][1] = pfac*(2*r[0]*r[1]);
        Dinsl[d][2] = pfac*(2*r[0]*r[2]);
        Dinsl[d][3] = pfac*(2*r[1]*r[1] + r[0]*r[0] - r2);
        Dinsl[d][4] = pfac*(2*r[1]*r[2]);
<<<<<<< HEAD
        
        if (ms)
        {
            for(i=0; i<5; i++)
=======

        if (ms)
        {
            for (i = 0; i < 5; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                Dins[d][i] = Dinsl[d][i]*invn;
            }
        }

        d++;
    }
<<<<<<< HEAD
  
    if (ms)
    {
        gmx_sum_sim(5*od->nr,Dins[0],ms);
    }
   
    /* Calculate the order tensor S for each experiment via optimization */
    for(ex=0; ex<od->nex; ex++)
    {
        for(i=0; i<5; i++)
        {
            rhs[ex][i] = 0;
            for(j=0; j<=i; j++)
=======

    if (ms)
    {
        gmx_sum_sim(5*od->nr, Dins[0], ms);
    }

    /* Calculate the order tensor S for each experiment via optimization */
    for (ex = 0; ex < od->nex; ex++)
    {
        for (i = 0; i < 5; i++)
        {
            rhs[ex][i] = 0;
            for (j = 0; j <= i; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                T[ex][i][j] = 0;
            }
        }
    }
    d = 0;
<<<<<<< HEAD
    for(fa=0; fa<nfa; fa+=3)
=======
    for (fa = 0; fa < nfa; fa += 3)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (bTAV)
        {
            /* Here we update Dtav in t_fcdata using the data in history_t.
             * Thus the results stay correct when this routine
             * is called multiple times.
             */
<<<<<<< HEAD
            for(i=0; i<5; i++)
=======
            for (i = 0; i < 5; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                Dtav[d][i] = edt*hist->orire_Dtav[d*5+i] + edt_1*Dins[d][i];
            }
        }
<<<<<<< HEAD
        
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        type   = forceatoms[fa];
        ex     = ip[type].orires.ex;
        weight = ip[type].orires.kfac;
        /* Calculate the vector rhs and half the matrix T for the 5 equations */
<<<<<<< HEAD
        for(i=0; i<5; i++) {
            rhs[ex][i] += Dtav[d][i]*ip[type].orires.obs*weight;
            for(j=0; j<=i; j++)
=======
        for (i = 0; i < 5; i++)
        {
            rhs[ex][i] += Dtav[d][i]*ip[type].orires.obs*weight;
            for (j = 0; j <= i; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                T[ex][i][j] += Dtav[d][i]*Dtav[d][j]*weight;
            }
        }
        d++;
    }
    /* Now we have all the data we can calculate S */
<<<<<<< HEAD
    for(ex=0; ex<od->nex; ex++)
    {
        /* Correct corrfac and copy one half of T to the other half */
        for(i=0; i<5; i++)
        {
            rhs[ex][i]  *= corrfac;
            T[ex][i][i] *= sqr(corrfac);
            for(j=0; j<i; j++)
=======
    for (ex = 0; ex < od->nex; ex++)
    {
        /* Correct corrfac and copy one half of T to the other half */
        for (i = 0; i < 5; i++)
        {
            rhs[ex][i]  *= corrfac;
            T[ex][i][i] *= sqr(corrfac);
            for (j = 0; j < i; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                T[ex][i][j] *= sqr(corrfac);
                T[ex][j][i]  = T[ex][i][j];
            }
        }
<<<<<<< HEAD
        m_inv_gen(T[ex],5,T[ex]);
=======
        m_inv_gen(T[ex], 5, T[ex]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Calculate the orientation tensor S for this experiment */
        S[ex][0][0] = 0;
        S[ex][0][1] = 0;
        S[ex][0][2] = 0;
        S[ex][1][1] = 0;
        S[ex][1][2] = 0;
<<<<<<< HEAD
        for(i=0; i<5; i++)
=======
        for (i = 0; i < 5; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            S[ex][0][0] += 1.5*T[ex][0][i]*rhs[ex][i];
            S[ex][0][1] += 1.5*T[ex][1][i]*rhs[ex][i];
            S[ex][0][2] += 1.5*T[ex][2][i]*rhs[ex][i];
            S[ex][1][1] += 1.5*T[ex][3][i]*rhs[ex][i];
            S[ex][1][2] += 1.5*T[ex][4][i]*rhs[ex][i];
        }
        S[ex][1][0] = S[ex][0][1];
        S[ex][2][0] = S[ex][0][2];
        S[ex][2][1] = S[ex][1][2];
        S[ex][2][2] = -S[ex][0][0] - S[ex][1][1];
    }
<<<<<<< HEAD
    
    wsv2 = 0;
    sw   = 0;
    
    d = 0;
    for(fa=0; fa<nfa; fa+=3)
    {
        type = forceatoms[fa];
        ex = ip[type].orires.ex;
        
=======

    wsv2 = 0;
    sw   = 0;

    d = 0;
    for (fa = 0; fa < nfa; fa += 3)
    {
        type = forceatoms[fa];
        ex   = ip[type].orires.ex;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        od->otav[d] = two_thr*
            corrfac*(S[ex][0][0]*Dtav[d][0] + S[ex][0][1]*Dtav[d][1] +
                     S[ex][0][2]*Dtav[d][2] + S[ex][1][1]*Dtav[d][3] +
                     S[ex][1][2]*Dtav[d][4]);
        if (bTAV)
        {
            od->oins[d] = two_thr*(S[ex][0][0]*Dins[d][0] + S[ex][0][1]*Dins[d][1] +
                                   S[ex][0][2]*Dins[d][2] + S[ex][1][1]*Dins[d][3] +
                                   S[ex][1][2]*Dins[d][4]);
        }
        if (ms)
        {
            /* When ensemble averaging is used recalculate the local orientation
             * for output to the energy file.
             */
            od->oinsl[d] = two_thr*
                (S[ex][0][0]*Dinsl[d][0] + S[ex][0][1]*Dinsl[d][1] +
                 S[ex][0][2]*Dinsl[d][2] + S[ex][1][1]*Dinsl[d][3] +
                 S[ex][1][2]*Dinsl[d][4]);
        }
<<<<<<< HEAD
        
        dev = od->otav[d] - ip[type].orires.obs;
        
        wsv2 += ip[type].orires.kfac*sqr(dev);
        sw   += ip[type].orires.kfac;
        
        d++;
    }
    od->rmsdev = sqrt(wsv2/sw);
    
    /* Rotate the S matrices back, so we get the correct grad(tr(S D)) */
    for(ex=0; ex<od->nex; ex++)
    {
        tmmul(R,S[ex],TMP);
        mmul(TMP,R,S[ex]);
    }
    
    return od->rmsdev;
    
    /* Approx. 120*nfa/3 flops */
}

real orires(int nfa,const t_iatom forceatoms[],const t_iparams ip[],
            const rvec x[],rvec f[],rvec fshift[],
            const t_pbc *pbc,const t_graph *g,
            real lambda,real *dvdlambda,
            const t_mdatoms *md,t_fcdata *fcd,
            int *global_atom_index)
{
    atom_id      ai,aj;
    int          fa,d,i,type,ex,power,ki=CENTRAL;
    ivec         dt;
    real         r2,invr,invr2,fc,smooth_fc,dev,devins,pfac;
    rvec         r,Sr,fij;
    real         vtot;
    const t_oriresdata *od;
    gmx_bool         bTAV;
    
    vtot = 0;
    od = &(fcd->orires);
    
=======

        dev = od->otav[d] - ip[type].orires.obs;

        wsv2 += ip[type].orires.kfac*sqr(dev);
        sw   += ip[type].orires.kfac;

        d++;
    }
    od->rmsdev = sqrt(wsv2/sw);

    /* Rotate the S matrices back, so we get the correct grad(tr(S D)) */
    for (ex = 0; ex < od->nex; ex++)
    {
        tmmul(R, S[ex], TMP);
        mmul(TMP, R, S[ex]);
    }

    return od->rmsdev;

    /* Approx. 120*nfa/3 flops */
}

real orires(int nfa, const t_iatom forceatoms[], const t_iparams ip[],
            const rvec x[], rvec f[], rvec fshift[],
            const t_pbc *pbc, const t_graph *g,
            real lambda, real *dvdlambda,
            const t_mdatoms *md, t_fcdata *fcd,
            int *global_atom_index)
{
    atom_id             ai, aj;
    int                 fa, d, i, type, ex, power, ki = CENTRAL;
    ivec                dt;
    real                r2, invr, invr2, fc, smooth_fc, dev, devins, pfac;
    rvec                r, Sr, fij;
    real                vtot;
    const t_oriresdata *od;
    gmx_bool            bTAV;

    vtot = 0;
    od   = &(fcd->orires);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (od->fc != 0)
    {
        bTAV = (od->edt != 0);

        smooth_fc = od->fc;
        if (bTAV)
        {
            /* Smoothly switch on the restraining when time averaging is used */
            smooth_fc *= (1.0 - od->exp_min_t_tau);
        }
<<<<<<< HEAD
        
        d = 0;
        for(fa=0; fa<nfa; fa+=3)
=======

        d = 0;
        for (fa = 0; fa < nfa; fa += 3)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            type  = forceatoms[fa];
            ai    = forceatoms[fa+1];
            aj    = forceatoms[fa+2];
            if (pbc)
            {
<<<<<<< HEAD
                ki = pbc_dx_aiuc(pbc,x[ai],x[aj],r);
            }
            else
            {
                rvec_sub(x[ai],x[aj],r);
=======
                ki = pbc_dx_aiuc(pbc, x[ai], x[aj], r);
            }
            else
            {
                rvec_sub(x[ai], x[aj], r);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            r2    = norm2(r);
            invr  = gmx_invsqrt(r2);
            invr2 = invr*invr;
            ex    = ip[type].orires.ex;
            power = ip[type].orires.power;
            fc    = smooth_fc*ip[type].orires.kfac;
            dev   = od->otav[d] - ip[type].orires.obs;
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* NOTE:
             * there is no real potential when time averaging is applied
             */
            vtot += 0.5*fc*sqr(dev);
<<<<<<< HEAD
            
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (bTAV)
            {
                /* Calculate the force as the sqrt of tav times instantaneous */
                devins = od->oins[d] - ip[type].orires.obs;
                if (dev*devins <= 0)
                {
                    dev = 0;
                }
                else
                {
                    dev = sqrt(dev*devins);
                    if (devins < 0)
                    {
                        dev = -dev;
                    }
                }
            }
<<<<<<< HEAD
            
            pfac  = fc*ip[type].orires.c*invr2;
            for(i=0; i<power; i++)
            {
                pfac *= invr;
            }
            mvmul(od->S[ex],r,Sr);
            for(i=0; i<DIM; i++)
            {
                fij[i] =
                    -pfac*dev*(4*Sr[i] - 2*(2+power)*invr2*iprod(Sr,r)*r[i]);
            }
            
            if (g)
            {
                ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
                ki=IVEC2IS(dt);
            }
            
            for(i=0; i<DIM; i++)
=======

            pfac  = fc*ip[type].orires.c*invr2;
            for (i = 0; i < power; i++)
            {
                pfac *= invr;
            }
            mvmul(od->S[ex], r, Sr);
            for (i = 0; i < DIM; i++)
            {
                fij[i] =
                    -pfac*dev*(4*Sr[i] - 2*(2+power)*invr2*iprod(Sr, r)*r[i]);
            }

            if (g)
            {
                ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
                ki = IVEC2IS(dt);
            }

            for (i = 0; i < DIM; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                f[ai][i]           += fij[i];
                f[aj][i]           -= fij[i];
                fshift[ki][i]      += fij[i];
                fshift[CENTRAL][i] -= fij[i];
            }
            d++;
        }
    }
<<<<<<< HEAD
    
    return vtot;
    
    /* Approx. 80*nfa/3 flops */
}

void update_orires_history(t_fcdata *fcd,history_t *hist)
{
    t_oriresdata *od;
    int pair,i;
    
=======

    return vtot;

    /* Approx. 80*nfa/3 flops */
}

void update_orires_history(t_fcdata *fcd, history_t *hist)
{
    t_oriresdata *od;
    int           pair, i;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    od = &(fcd->orires);
    if (od->edt != 0)
    {
        /* Copy the new time averages that have been calculated
         *  in calc_orires_dev.
         */
        hist->orire_initf = od->exp_min_t_tau;
<<<<<<< HEAD
        for(pair=0; pair<od->nr; pair++)
        {
            for(i=0; i<5; i++)
=======
        for (pair = 0; pair < od->nr; pair++)
        {
            for (i = 0; i < 5; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                hist->orire_Dtav[pair*5+i] = od->Dtav[pair][i];
            }
        }
    }
}
