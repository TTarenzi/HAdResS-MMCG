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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

<<<<<<< HEAD
#ifdef GMX_THREAD_SHM_FDECOMP
#include <pthread.h> 
#endif

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "network.h"
#include "nsgrid.h"
#include "force.h"
#include "nonbonded.h"
#include "ns.h"
#include "pbc.h"
#include "names.h"
#include "gmx_fatal.h"
#include "nrnb.h"
#include "txtdump.h"
#include "mtop_util.h"

#include "domdec.h"
#include "adress.h"


<<<<<<< HEAD
/* 
=======
/*
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *    E X C L U S I O N   H A N D L I N G
 */

#ifdef DEBUG
<<<<<<< HEAD
static void SETEXCL_(t_excl e[],atom_id i,atom_id j)
{   e[j] = e[j] | (1<<i); }
static void RMEXCL_(t_excl e[],atom_id i,atom_id j) 
{ e[j]=e[j] & ~(1<<i); }
static gmx_bool ISEXCL_(t_excl e[],atom_id i,atom_id j) 
{ return (gmx_bool)(e[j] & (1<<i)); }
static gmx_bool NOTEXCL_(t_excl e[],atom_id i,atom_id j)
{  return !(ISEXCL(e,i,j)); }
#else
#define SETEXCL(e,i,j) (e)[((atom_id) (j))] |= (1<<((atom_id) (i)))
#define RMEXCL(e,i,j)  (e)[((atom_id) (j))] &= (~(1<<((atom_id) (i))))
#define ISEXCL(e,i,j)  (gmx_bool) ((e)[((atom_id) (j))] & (1<<((atom_id) (i))))
#define NOTEXCL(e,i,j) !(ISEXCL(e,i,j))
#endif

=======
static void SETEXCL_(t_excl e[], atom_id i, atom_id j)
{
    e[j] = e[j] | (1<<i);
}
static void RMEXCL_(t_excl e[], atom_id i, atom_id j)
{
    e[j] = e[j] & ~(1<<i);
}
static gmx_bool ISEXCL_(t_excl e[], atom_id i, atom_id j)
{
    return (gmx_bool)(e[j] & (1<<i));
}
static gmx_bool NOTEXCL_(t_excl e[], atom_id i, atom_id j)
{
    return !(ISEXCL(e, i, j));
}
#else
#define SETEXCL(e, i, j) (e)[((atom_id) (j))] |= (1<<((atom_id) (i)))
#define RMEXCL(e, i, j)  (e)[((atom_id) (j))] &= (~(1<<((atom_id) (i))))
#define ISEXCL(e, i, j)  (gmx_bool) ((e)[((atom_id) (j))] & (1<<((atom_id) (i))))
#define NOTEXCL(e, i, j) !(ISEXCL(e, i, j))
#endif

static int
round_up_to_simd_width(int length, int simd_width)
{
    int offset, newlength;

    offset = (simd_width > 0) ? length % simd_width : 0;

    return (offset == 0) ? length : length-offset+simd_width;
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/************************************************
 *
 *  U T I L I T I E S    F O R    N S
 *
 ************************************************/

static void reallocate_nblist(t_nblist *nl)
{
    if (gmx_debug_at)
    {
<<<<<<< HEAD
        fprintf(debug,"reallocating neigborlist il_code=%d, maxnri=%d\n",
                nl->il_code,nl->maxnri); 
    }
    srenew(nl->iinr,   nl->maxnri);
    if (nl->enlist == enlistCG_CG)
    {
        srenew(nl->iinr_end,nl->maxnri);
=======
        fprintf(debug, "reallocating neigborlist (ielec=%d, ivdw=%d, igeometry=%d, type=%d), maxnri=%d\n",
                nl->ielec, nl->ivdw, nl->igeometry, nl->type, nl->maxnri);
    }
    srenew(nl->iinr,   nl->maxnri);
    if (nl->igeometry == GMX_NBLIST_GEOMETRY_CG_CG)
    {
        srenew(nl->iinr_end, nl->maxnri);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    srenew(nl->gid,    nl->maxnri);
    srenew(nl->shift,  nl->maxnri);
    srenew(nl->jindex, nl->maxnri+1);
}

<<<<<<< HEAD
/* ivdw/icoul are used to determine the type of interaction, so we
 * can set an innerloop index here. The obvious choice for this would have
 * been the vdwtype/coultype values in the forcerecord, but unfortunately 
 * those types are braindead - for instance both Buckingham and normal 
 * Lennard-Jones use the same value (evdwCUT), and a separate gmx_boolean variable
 * to determine which interaction is used. There is further no special value
 * for 'no interaction'. For backward compatibility with old TPR files we won't
 * change this in the 3.x series, so when calling this routine you should use:
 *
 * icoul=0 no coulomb interaction
 * icoul=1 cutoff standard coulomb
 * icoul=2 reaction-field coulomb
 * icoul=3 tabulated coulomb
 *
 * ivdw=0 no vdw interaction
 * ivdw=1 standard L-J interaction
 * ivdw=2 Buckingham
 * ivdw=3 tabulated vdw.
 *
 * Kind of ugly, but it works.
 */
static void init_nblist(t_nblist *nl_sr,t_nblist *nl_lr,
                        int maxsr,int maxlr,
                        int ivdw, int icoul, 
                        gmx_bool bfree, int enlist)
{
    t_nblist *nl;
    int      homenr;
    int      i,nn;
    
    int inloop[20] =
    { 
        eNR_NBKERNEL_NONE,
        eNR_NBKERNEL010,
        eNR_NBKERNEL020,
        eNR_NBKERNEL030,
        eNR_NBKERNEL100,
        eNR_NBKERNEL110,
        eNR_NBKERNEL120,
        eNR_NBKERNEL130,
        eNR_NBKERNEL200,
        eNR_NBKERNEL210,
        eNR_NBKERNEL220,
        eNR_NBKERNEL230,
        eNR_NBKERNEL300,
        eNR_NBKERNEL310,
        eNR_NBKERNEL320,
        eNR_NBKERNEL330,
        eNR_NBKERNEL400,
        eNR_NBKERNEL410,
        eNR_NBKERNEL_NONE,
        eNR_NBKERNEL430
    };
  
    for(i=0; (i<2); i++)
=======

static void init_nblist(FILE *log, t_nblist *nl_sr, t_nblist *nl_lr,
                        int maxsr, int maxlr,
                        int ivdw, int ivdwmod,
                        int ielec, int ielecmod,
                        int igeometry, int type)
{
    t_nblist *nl;
    int       homenr;
    int       i, nn;

    for (i = 0; (i < 2); i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        nl     = (i == 0) ? nl_sr : nl_lr;
        homenr = (i == 0) ? maxsr : maxlr;

        if (nl == NULL)
        {
            continue;
        }
<<<<<<< HEAD
        
        /* Set coul/vdw in neighborlist, and for the normal loops we determine
         * an index of which one to call.
         */
        nl->ivdw  = ivdw;
        nl->icoul = icoul;
        nl->free_energy = bfree;
    
        if (bfree)
        {
            nl->enlist  = enlistATOM_ATOM;
            nl->il_code = eNR_NBKERNEL_FREE_ENERGY;
        }
        else
        {
            nl->enlist = enlist;

            nn = inloop[4*icoul + ivdw];
            
            /* solvent loops follow directly after the corresponding
            * ordinary loops, in the order:
            *
            * SPC, SPC-SPC, TIP4p, TIP4p-TIP4p
            *   
            */
            switch (enlist) {
            case enlistATOM_ATOM:
            case enlistCG_CG:
                break;
            case enlistSPC_ATOM:     nn += 1; break;
            case enlistSPC_SPC:      nn += 2; break;
            case enlistTIP4P_ATOM:   nn += 3; break;
            case enlistTIP4P_TIP4P:  nn += 4; break;
            }
            
            nl->il_code = nn;
        }

        if (debug)
            fprintf(debug,"Initiating neighbourlist type %d for %s interactions,\nwith %d SR, %d LR atoms.\n",
                    nl->il_code,ENLISTTYPE(enlist),maxsr,maxlr);
        
=======


        /* Set coul/vdw in neighborlist, and for the normal loops we determine
         * an index of which one to call.
         */
        nl->ivdw        = ivdw;
        nl->ivdwmod     = ivdwmod;
        nl->ielec       = ielec;
        nl->ielecmod    = ielecmod;
        nl->type        = type;
        nl->igeometry   = igeometry;

        if (nl->type == GMX_NBLIST_INTERACTION_FREE_ENERGY)
        {
            nl->igeometry  = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
        }

        /* This will also set the simd_padding_width field */
        gmx_nonbonded_set_kernel_pointers( (i == 0) ? log : NULL, nl);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* maxnri is influenced by the number of shifts (maximum is 8)
         * and the number of energy groups.
         * If it is not enough, nl memory will be reallocated during the run.
         * 4 seems to be a reasonable factor, which only causes reallocation
         * during runs with tiny and many energygroups.
         */
        nl->maxnri      = homenr*4;
        nl->maxnrj      = 0;
        nl->maxlen      = 0;
        nl->nri         = -1;
        nl->nrj         = 0;
        nl->iinr        = NULL;
        nl->gid         = NULL;
        nl->shift       = NULL;
        nl->jindex      = NULL;
        reallocate_nblist(nl);
        nl->jindex[0] = 0;
<<<<<<< HEAD
#ifdef GMX_THREAD_SHM_FDECOMP
        nl->counter = 0;
        snew(nl->mtx,1);
        pthread_mutex_init(nl->mtx,NULL);
#endif
    }
}

void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr)
{
   /* Make maxlr tunable! (does not seem to be a big difference though) 
    * This parameter determines the number of i particles in a long range 
    * neighbourlist. Too few means many function calls, too many means
    * cache trashing.
    */
   int maxsr,maxsr_wat,maxlr,maxlr_wat;
   int icoul,icoulf,ivdw;
   int solvent;
   int enlist_def,enlist_w,enlist_ww;
   int i;
   t_nblists *nbl;

   /* maxsr     = homenr-fr->nWatMol*3; */
   maxsr     = homenr;

   if (maxsr < 0)
   {
     gmx_fatal(FARGS,"%s, %d: Negative number of short range atoms.\n"
		 "Call your Gromacs dealer for assistance.",__FILE__,__LINE__);
   }
   /* This is just for initial allocation, so we do not reallocate
    * all the nlist arrays many times in a row.
    * The numbers seem very accurate, but they are uncritical.
    */
   maxsr_wat = min(fr->nWatMol,(homenr+2)/3); 
   if (fr->bTwinRange) 
   {
       maxlr     = 50;
       maxlr_wat = min(maxsr_wat,maxlr);
   }
   else
   {
     maxlr = maxlr_wat = 0;
   }  

   /* Determine the values for icoul/ivdw. */
   /* Start with GB */
   if(fr->bGB)
   {
       icoul=4;
   }
   else if (fr->bcoultab)
   {
       icoul = 3;
   }
   else if (EEL_RF(fr->eeltype))
   {
       icoul = 2;
   }
   else 
   {
       icoul = 1;
   }
   
   if (fr->bvdwtab)
   {
       ivdw = 3;
   }
   else if (fr->bBHAM)
   {
       ivdw = 2;
   }
   else 
   {
       ivdw = 1;
   }

   fr->ns.bCGlist = (getenv("GMX_NBLISTCG") != 0);
   if (!fr->ns.bCGlist)
   {
       enlist_def = enlistATOM_ATOM;
   }
   else
   {
       enlist_def = enlistCG_CG;
       if (log != NULL)
       {
           fprintf(log,"\nUsing charge-group - charge-group neighbor lists and kernels\n\n");
       }
   }
   
   if (fr->solvent_opt == esolTIP4P) {
       enlist_w  = enlistTIP4P_ATOM;
       enlist_ww = enlistTIP4P_TIP4P;
   } else {
       enlist_w  = enlistSPC_ATOM;
       enlist_ww = enlistSPC_SPC;
   }

   for(i=0; i<fr->nnblists; i++) 
   {
       nbl = &(fr->nblists[i]);
       init_nblist(&nbl->nlist_sr[eNL_VDWQQ],&nbl->nlist_lr[eNL_VDWQQ],
                   maxsr,maxlr,ivdw,icoul,FALSE,enlist_def);
       init_nblist(&nbl->nlist_sr[eNL_VDW],&nbl->nlist_lr[eNL_VDW],
                   maxsr,maxlr,ivdw,0,FALSE,enlist_def);
       init_nblist(&nbl->nlist_sr[eNL_QQ],&nbl->nlist_lr[eNL_QQ],
                   maxsr,maxlr,0,icoul,FALSE,enlist_def);
       init_nblist(&nbl->nlist_sr[eNL_VDWQQ_WATER],&nbl->nlist_lr[eNL_VDWQQ_WATER],
                   maxsr_wat,maxlr_wat,ivdw,icoul, FALSE,enlist_w);
       init_nblist(&nbl->nlist_sr[eNL_QQ_WATER],&nbl->nlist_lr[eNL_QQ_WATER],
                   maxsr_wat,maxlr_wat,0,icoul, FALSE,enlist_w);
       init_nblist(&nbl->nlist_sr[eNL_VDWQQ_WATERWATER],&nbl->nlist_lr[eNL_VDWQQ_WATERWATER],
                   maxsr_wat,maxlr_wat,ivdw,icoul, FALSE,enlist_ww);
       init_nblist(&nbl->nlist_sr[eNL_QQ_WATERWATER],&nbl->nlist_lr[eNL_QQ_WATERWATER],
                   maxsr_wat,maxlr_wat,0,icoul, FALSE,enlist_ww);
       
       if (fr->efep != efepNO) 
       {
           if (fr->bEwald)
           {
               icoulf = 5;
           }
           else
           {
               icoulf = icoul;
           }

           init_nblist(&nbl->nlist_sr[eNL_VDWQQ_FREE],&nbl->nlist_lr[eNL_VDWQQ_FREE],
                       maxsr,maxlr,ivdw,icoulf,TRUE,enlistATOM_ATOM);
           init_nblist(&nbl->nlist_sr[eNL_VDW_FREE],&nbl->nlist_lr[eNL_VDW_FREE],
                       maxsr,maxlr,ivdw,0,TRUE,enlistATOM_ATOM);
           init_nblist(&nbl->nlist_sr[eNL_QQ_FREE],&nbl->nlist_lr[eNL_QQ_FREE],
                       maxsr,maxlr,0,icoulf,TRUE,enlistATOM_ATOM);
       }  
   }
   /* QMMM MM list */
   if (fr->bQMMM && fr->qr->QMMMscheme != eQMMMschemeoniom)
   {
       init_nblist(&fr->QMMMlist,NULL,
                   maxsr,maxlr,0,icoul,FALSE,enlistATOM_ATOM);
   }

   fr->ns.nblist_initialized=TRUE;
=======

        if (debug)
        {
            fprintf(debug, "Initiating neighbourlist (ielec=%d, ivdw=%d, type=%d) for %s interactions,\nwith %d SR, %d LR atoms.\n",
                    nl->ielec, nl->ivdw, nl->type, gmx_nblist_geometry_names[nl->igeometry], maxsr, maxlr);
        }
    }
}

void init_neighbor_list(FILE *log, t_forcerec *fr, int homenr)
{
    /* Make maxlr tunable! (does not seem to be a big difference though)
     * This parameter determines the number of i particles in a long range
     * neighbourlist. Too few means many function calls, too many means
     * cache trashing.
     */
    int        maxsr, maxsr_wat, maxlr, maxlr_wat;
    int        ielec, ielecf, ivdw, ielecmod, ielecmodf, ivdwmod, type;
    int        solvent;
    int        igeometry_def, igeometry_w, igeometry_ww;
    int        i;
    t_nblists *nbl;

    /* maxsr     = homenr-fr->nWatMol*3; */
    maxsr     = homenr;

    if (maxsr < 0)
    {
        gmx_fatal(FARGS, "%s, %d: Negative number of short range atoms.\n"
                  "Call your Gromacs dealer for assistance.", __FILE__, __LINE__);
    }
    /* This is just for initial allocation, so we do not reallocate
     * all the nlist arrays many times in a row.
     * The numbers seem very accurate, but they are uncritical.
     */
    maxsr_wat = min(fr->nWatMol, (homenr+2)/3);
    if (fr->bTwinRange)
    {
        maxlr     = 50;
        maxlr_wat = min(maxsr_wat, maxlr);
    }
    else
    {
        maxlr = maxlr_wat = 0;
    }

    /* Determine the values for ielec/ivdw. */
    ielec    = fr->nbkernel_elec_interaction;
    ivdw     = fr->nbkernel_vdw_interaction;
    ielecmod = fr->nbkernel_elec_modifier;
    ivdwmod  = fr->nbkernel_vdw_modifier;
    type     = GMX_NBLIST_INTERACTION_STANDARD;

    fr->ns.bCGlist = (getenv("GMX_NBLISTCG") != 0);
    if (!fr->ns.bCGlist)
    {
        igeometry_def = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
    }
    else
    {
        igeometry_def = GMX_NBLIST_GEOMETRY_CG_CG;
        if (log != NULL)
        {
            fprintf(log, "\nUsing charge-group - charge-group neighbor lists and kernels\n\n");
        }
    }

    if (fr->solvent_opt == esolTIP4P)
    {
        igeometry_w  = GMX_NBLIST_GEOMETRY_WATER4_PARTICLE;
        igeometry_ww = GMX_NBLIST_GEOMETRY_WATER4_WATER4;
    }
    else
    {
        igeometry_w  = GMX_NBLIST_GEOMETRY_WATER3_PARTICLE;
        igeometry_ww = GMX_NBLIST_GEOMETRY_WATER3_WATER3;
    }

    for (i = 0; i < fr->nnblists; i++)
    {
        nbl = &(fr->nblists[i]);

        if ((fr->adress_type != eAdressOff) && (i >= fr->nnblists/2))
        {
            type = GMX_NBLIST_INTERACTION_ADRESS;
        }
        init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ], &nbl->nlist_lr[eNL_VDWQQ],
                    maxsr, maxlr, ivdw, ivdwmod, ielec, ielecmod, igeometry_def, type);
        init_nblist(log, &nbl->nlist_sr[eNL_VDW], &nbl->nlist_lr[eNL_VDW],
                    maxsr, maxlr, ivdw, ivdwmod, GMX_NBKERNEL_ELEC_NONE, eintmodNONE, igeometry_def, type);
        init_nblist(log, &nbl->nlist_sr[eNL_QQ], &nbl->nlist_lr[eNL_QQ],
                    maxsr, maxlr, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, igeometry_def, type);
        init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ_WATER], &nbl->nlist_lr[eNL_VDWQQ_WATER],
                    maxsr_wat, maxlr_wat, ivdw, ivdwmod, ielec, ielecmod, igeometry_w, type);
        init_nblist(log, &nbl->nlist_sr[eNL_QQ_WATER], &nbl->nlist_lr[eNL_QQ_WATER],
                    maxsr_wat, maxlr_wat, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, igeometry_w, type);
        init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ_WATERWATER], &nbl->nlist_lr[eNL_VDWQQ_WATERWATER],
                    maxsr_wat, maxlr_wat, ivdw, ivdwmod, ielec, ielecmod, igeometry_ww, type);
        init_nblist(log, &nbl->nlist_sr[eNL_QQ_WATERWATER], &nbl->nlist_lr[eNL_QQ_WATERWATER],
                    maxsr_wat, maxlr_wat, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, igeometry_ww, type);

        /* Did we get the solvent loops so we can use optimized water kernels? */
        if (nbl->nlist_sr[eNL_VDWQQ_WATER].kernelptr_vf == NULL
            || nbl->nlist_sr[eNL_QQ_WATER].kernelptr_vf == NULL
#ifndef DISABLE_WATERWATER_NLIST
            || nbl->nlist_sr[eNL_VDWQQ_WATERWATER].kernelptr_vf == NULL
            || nbl->nlist_sr[eNL_QQ_WATERWATER].kernelptr_vf == NULL
#endif
            )
        {
            fr->solvent_opt = esolNO;
            fprintf(log, "Note: The available nonbonded kernels do not support water optimization - disabling.\n");
        }

        if (fr->efep != efepNO)
        {
            if ((fr->bEwald) && (fr->sc_alphacoul > 0)) /* need to handle long range differently if using softcore */
            {
                ielecf    = GMX_NBKERNEL_ELEC_EWALD;
                ielecmodf = eintmodNONE;
            }
            else
            {
                ielecf    = ielec;
                ielecmodf = ielecmod;
            }

            init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ_FREE], &nbl->nlist_lr[eNL_VDWQQ_FREE],
                        maxsr, maxlr, ivdw, ivdwmod, ielecf, ielecmod, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_FREE_ENERGY);
            init_nblist(log, &nbl->nlist_sr[eNL_VDW_FREE], &nbl->nlist_lr[eNL_VDW_FREE],
                        maxsr, maxlr, ivdw, ivdwmod, GMX_NBKERNEL_ELEC_NONE, eintmodNONE, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_FREE_ENERGY);
            init_nblist(log, &nbl->nlist_sr[eNL_QQ_FREE], &nbl->nlist_lr[eNL_QQ_FREE],
                        maxsr, maxlr, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielecf, ielecmod, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_FREE_ENERGY);
        }
    }
    /* QMMM MM list */
    if (fr->bQMMM && fr->qr->QMMMscheme != eQMMMschemeoniom)
    {
        init_nblist(log, &fr->QMMMlist, NULL,
                    maxsr, maxlr, 0, 0, ielec, ielecmod, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_STANDARD);
    }

    if (log != NULL)
    {
        fprintf(log, "\n");
    }

    fr->ns.nblist_initialized = TRUE;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

static void reset_nblist(t_nblist *nl)
{
<<<<<<< HEAD
     nl->nri       = -1;
     nl->nrj       = 0;
     nl->maxlen    = 0;
     if (nl->jindex)
     {
         nl->jindex[0] = 0;
     }
}

static void reset_neighbor_list(t_forcerec *fr,gmx_bool bLR,int nls,int eNL)
{
    int n,i;
  
    if (bLR) 
    {
        reset_nblist(&(fr->nblists[nls].nlist_lr[eNL]));
    }
    else 
    {
        for(n=0; n<fr->nnblists; n++)
        {
            for(i=0; i<eNL_NR; i++)
            {
                reset_nblist(&(fr->nblists[n].nlist_sr[i]));
            }
        }
        if (fr->bQMMM)
        { 
            /* only reset the short-range nblist */
            reset_nblist(&(fr->QMMMlist));
=======
    nl->nri       = -1;
    nl->nrj       = 0;
    nl->maxlen    = 0;
    if (nl->jindex)
    {
        nl->jindex[0] = 0;
    }
}

static void reset_neighbor_lists(t_forcerec *fr, gmx_bool bResetSR, gmx_bool bResetLR)
{
    int n, i;

    if (fr->bQMMM)
    {
        /* only reset the short-range nblist */
        reset_nblist(&(fr->QMMMlist));
    }

    for (n = 0; n < fr->nnblists; n++)
    {
        for (i = 0; i < eNL_NR; i++)
        {
            if (bResetSR)
            {
                reset_nblist( &(fr->nblists[n].nlist_sr[i]) );
            }
            if (bResetLR)
            {
                reset_nblist( &(fr->nblists[n].nlist_lr[i]) );
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
}




static inline void new_i_nblist(t_nblist *nlist,
<<<<<<< HEAD
                                gmx_bool bLR,atom_id i_atom,int shift,int gid)
{
    int    i,k,nri,nshift;
    
    nri = nlist->nri;
    
    /* Check whether we have to increase the i counter */
    if ((nri == -1) ||
        (nlist->iinr[nri]  != i_atom) || 
        (nlist->shift[nri] != shift) || 
        (nlist->gid[nri]   != gid))
    {
        /* This is something else. Now see if any entries have 
         * been added in the list of the previous atom.
         */
        if ((nri == -1) ||
            ((nlist->jindex[nri+1] > nlist->jindex[nri]) && 
=======
                                gmx_bool bLR, atom_id i_atom, int shift, int gid)
{
    int    i, k, nri, nshift;

    nri = nlist->nri;

    /* Check whether we have to increase the i counter */
    if ((nri == -1) ||
        (nlist->iinr[nri]  != i_atom) ||
        (nlist->shift[nri] != shift) ||
        (nlist->gid[nri]   != gid))
    {
        /* This is something else. Now see if any entries have
         * been added in the list of the previous atom.
         */
        if ((nri == -1) ||
            ((nlist->jindex[nri+1] > nlist->jindex[nri]) &&
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
             (nlist->gid[nri] != -1)))
        {
            /* If so increase the counter */
            nlist->nri++;
            nri++;
            if (nlist->nri >= nlist->maxnri)
            {
                nlist->maxnri += over_alloc_large(nlist->nri);
                reallocate_nblist(nlist);
            }
        }
        /* Set the number of neighbours and the atom number */
        nlist->jindex[nri+1] = nlist->jindex[nri];
        nlist->iinr[nri]     = i_atom;
        nlist->gid[nri]      = gid;
        nlist->shift[nri]    = shift;
    }
}

<<<<<<< HEAD
static inline void close_i_nblist(t_nblist *nlist) 
{
    int nri = nlist->nri;
    int len;
    
    if (nri >= 0)
    {
        nlist->jindex[nri+1] = nlist->nrj;
        
        len=nlist->nrj -  nlist->jindex[nri];
        
        /* nlist length for water i molecules is treated statically 
         * in the innerloops 
=======
static inline void close_i_nblist(t_nblist *nlist)
{
    int nri = nlist->nri;
    int len;

    if (nri >= 0)
    {
        /* Add elements up to padding. Since we allocate memory in units
         * of the simd_padding width, we do not have to check for possible
         * list reallocation here.
         */
        while ((nlist->nrj % nlist->simd_padding_width) != 0)
        {
            /* Use -4 here, so we can write forces for 4 atoms before real data */
            nlist->jjnr[nlist->nrj++] = -4;
        }
        nlist->jindex[nri+1] = nlist->nrj;

        len = nlist->nrj -  nlist->jindex[nri];

        /* nlist length for water i molecules is treated statically
         * in the innerloops
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
         */
        if (len > nlist->maxlen)
        {
            nlist->maxlen = len;
        }
    }
}

static inline void close_nblist(t_nblist *nlist)
{
    /* Only close this nblist when it has been initialized.
     * Avoid the creation of i-lists with no j-particles.
     */
    if (nlist->nrj == 0)
    {
        /* Some assembly kernels do not support empty lists,
         * make sure here that we don't generate any empty lists.
         * With the current ns code this branch is taken in two cases:
         * No i-particles at all: nri=-1 here
         * There are i-particles, but no j-particles; nri=0 here
         */
        nlist->nri = 0;
    }
    else
    {
        /* Close list number nri by incrementing the count */
        nlist->nri++;
    }
}

<<<<<<< HEAD
static inline void close_neighbor_list(t_forcerec *fr,gmx_bool bLR,int nls,int eNL, 
                                       gmx_bool bMakeQMMMnblist)
{
    int n,i;
    
    if (bMakeQMMMnblist) {
        if (!bLR)
        {
            close_nblist(&(fr->QMMMlist));
        }
    }
    else 
    {
        if (bLR)
        {
            close_nblist(&(fr->nblists[nls].nlist_lr[eNL]));
        }
        else
        { 
            for(n=0; n<fr->nnblists; n++)
            {
                for(i=0; (i<eNL_NR); i++)
                {
                    close_nblist(&(fr->nblists[n].nlist_sr[i]));
                }
            }
=======
static inline void close_neighbor_lists(t_forcerec *fr, gmx_bool bMakeQMMMnblist)
{
    int n, i;

    if (bMakeQMMMnblist)
    {
        close_nblist(&(fr->QMMMlist));
    }

    for (n = 0; n < fr->nnblists; n++)
    {
        for (i = 0; (i < eNL_NR); i++)
        {
            close_nblist(&(fr->nblists[n].nlist_sr[i]));
            close_nblist(&(fr->nblists[n].nlist_lr[i]));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
}

<<<<<<< HEAD
static inline void add_j_to_nblist(t_nblist *nlist,atom_id j_atom,gmx_bool bLR)
{
    int nrj=nlist->nrj;
    
    if (nlist->nrj >= nlist->maxnrj)
    {
        nlist->maxnrj = over_alloc_small(nlist->nrj + 1);
        if (gmx_debug_at)
            fprintf(debug,"Increasing %s nblist %s j size to %d\n",
                    bLR ? "LR" : "SR",nrnb_str(nlist->il_code),nlist->maxnrj);
        
        srenew(nlist->jjnr,nlist->maxnrj);
    }

    nlist->jjnr[nrj] = j_atom;
    nlist->nrj ++;
}

static inline void add_j_to_nblist_cg(t_nblist *nlist,
                                      atom_id j_start,int j_end,
                                      t_excl *bexcl,gmx_bool i_is_j,
                                      gmx_bool bLR)
{
    int nrj=nlist->nrj;
=======

static inline void add_j_to_nblist(t_nblist *nlist, atom_id j_atom, gmx_bool bLR)
{
    int nrj = nlist->nrj;

    if (nlist->nrj >= nlist->maxnrj)
    {
        nlist->maxnrj = round_up_to_simd_width(over_alloc_small(nlist->nrj + 1), nlist->simd_padding_width);

        if (gmx_debug_at)
        {
            fprintf(debug, "Increasing %s nblist (ielec=%d,ivdw=%d,type=%d,igeometry=%d) j size to %d\n",
                    bLR ? "LR" : "SR", nlist->ielec, nlist->ivdw, nlist->type, nlist->igeometry, nlist->maxnrj);
        }

        srenew(nlist->jjnr, nlist->maxnrj);
    }

    nlist->jjnr[nrj] = j_atom;
    nlist->nrj++;
}

static inline void add_j_to_nblist_cg(t_nblist *nlist,
                                      atom_id j_start, int j_end,
                                      t_excl *bexcl, gmx_bool i_is_j,
                                      gmx_bool bLR)
{
    int nrj = nlist->nrj;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    int j;

    if (nlist->nrj >= nlist->maxnrj)
    {
        nlist->maxnrj = over_alloc_small(nlist->nrj + 1);
        if (gmx_debug_at)
<<<<<<< HEAD
            fprintf(debug,"Increasing %s nblist %s j size to %d\n",
                    bLR ? "LR" : "SR",nrnb_str(nlist->il_code),nlist->maxnrj);
        
        srenew(nlist->jjnr    ,nlist->maxnrj);
        srenew(nlist->jjnr_end,nlist->maxnrj);
        srenew(nlist->excl    ,nlist->maxnrj*MAX_CGCGSIZE);
=======
        {
            fprintf(debug, "Increasing %s nblist (ielec=%d,ivdw=%d,type=%d,igeometry=%d) j size to %d\n",
                    bLR ? "LR" : "SR", nlist->ielec, nlist->ivdw, nlist->type, nlist->igeometry, nlist->maxnrj);
        }

        srenew(nlist->jjnr, nlist->maxnrj);
        srenew(nlist->jjnr_end, nlist->maxnrj);
        srenew(nlist->excl, nlist->maxnrj*MAX_CGCGSIZE);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    nlist->jjnr[nrj]     = j_start;
    nlist->jjnr_end[nrj] = j_end;

    if (j_end - j_start > MAX_CGCGSIZE)
    {
<<<<<<< HEAD
        gmx_fatal(FARGS,"The charge-group - charge-group neighborlist do not support charge groups larger than %d, found a charge group of size %d",MAX_CGCGSIZE,j_end-j_start);
    }

    /* Set the exclusions */
    for(j=j_start; j<j_end; j++)
=======
        gmx_fatal(FARGS, "The charge-group - charge-group neighborlist do not support charge groups larger than %d, found a charge group of size %d", MAX_CGCGSIZE, j_end-j_start);
    }

    /* Set the exclusions */
    for (j = j_start; j < j_end; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        nlist->excl[nrj*MAX_CGCGSIZE + j - j_start] = bexcl[j];
    }
    if (i_is_j)
    {
        /* Avoid double counting of intra-cg interactions */
<<<<<<< HEAD
        for(j=1; j<j_end-j_start; j++)
=======
        for (j = 1; j < j_end-j_start; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            nlist->excl[nrj*MAX_CGCGSIZE + j] |= (1<<j) - 1;
        }
    }

<<<<<<< HEAD
    nlist->nrj ++;
}

typedef void
put_in_list_t(gmx_bool              bHaveVdW[],
              int               ngid,
              t_mdatoms *       md,
              int               icg,
              int               jgid,
              int               nj,
              atom_id           jjcg[],
              atom_id           index[],
              t_excl            bExcl[],
              int               shift,
              t_forcerec *      fr,
              gmx_bool              bLR,
              gmx_bool              bDoVdW,
              gmx_bool              bDoCoul);

static void 
put_in_list_at(gmx_bool              bHaveVdW[],
               int               ngid,
               t_mdatoms *       md,
               int               icg,
               int               jgid,
               int               nj,
               atom_id           jjcg[],
               atom_id           index[],
               t_excl            bExcl[],
               int               shift,
               t_forcerec *      fr,
               gmx_bool              bLR,
               gmx_bool              bDoVdW,
               gmx_bool              bDoCoul)
=======
    nlist->nrj++;
}

typedef void
    put_in_list_t (gmx_bool              bHaveVdW[],
                   int                   ngid,
                   t_mdatoms     *       md,
                   int                   icg,
                   int                   jgid,
                   int                   nj,
                   atom_id               jjcg[],
                   atom_id               index[],
                   t_excl                bExcl[],
                   int                   shift,
                   t_forcerec     *      fr,
                   gmx_bool              bLR,
                   gmx_bool              bDoVdW,
                   gmx_bool              bDoCoul,
                   int                   solvent_opt);

static void
put_in_list_at(gmx_bool              bHaveVdW[],
               int                   ngid,
               t_mdatoms     *       md,
               int                   icg,
               int                   jgid,
               int                   nj,
               atom_id               jjcg[],
               atom_id               index[],
               t_excl                bExcl[],
               int                   shift,
               t_forcerec     *      fr,
               gmx_bool              bLR,
               gmx_bool              bDoVdW,
               gmx_bool              bDoCoul,
               int                   solvent_opt)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    /* The a[] index has been removed,
     * to put it back in i_atom should be a[i0] and jj should be a[jj].
     */
<<<<<<< HEAD
    t_nblist *   vdwc;
    t_nblist *   vdw;
    t_nblist *   coul;
    t_nblist *   vdwc_free  = NULL;
    t_nblist *   vdw_free   = NULL;
    t_nblist *   coul_free  = NULL;
    t_nblist *   vdwc_ww    = NULL;
    t_nblist *   coul_ww    = NULL;
    
    int 	    i,j,jcg,igid,gid,nbl_ind,ind_ij;
    atom_id   jj,jj0,jj1,i_atom;
    int       i0,nicg,len;
    
    int       *cginfo;
    int       *type,*typeB;
    real      *charge,*chargeB;
    real      qi,qiB,qq,rlj;
    gmx_bool      bFreeEnergy,bFree,bFreeJ,bNotEx,*bPert;
    gmx_bool      bDoVdW_i,bDoCoul_i,bDoCoul_i_sol;
    int       iwater,jwater;
    t_nblist  *nlist;
    
=======
    t_nblist  *   vdwc;
    t_nblist  *   vdw;
    t_nblist  *   coul;
    t_nblist  *   vdwc_free  = NULL;
    t_nblist  *   vdw_free   = NULL;
    t_nblist  *   coul_free  = NULL;
    t_nblist  *   vdwc_ww    = NULL;
    t_nblist  *   coul_ww    = NULL;

    int           i, j, jcg, igid, gid, nbl_ind, ind_ij;
    atom_id       jj, jj0, jj1, i_atom;
    int           i0, nicg, len;

    int          *cginfo;
    int          *type, *typeB;
    real         *charge, *chargeB;
    real          qi, qiB, qq, rlj;
    gmx_bool      bFreeEnergy, bFree, bFreeJ, bNotEx, *bPert;
    gmx_bool      bDoVdW_i, bDoCoul_i, bDoCoul_i_sol;
    int           iwater, jwater;
    t_nblist     *nlist;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Copy some pointers */
    cginfo  = fr->cginfo;
    charge  = md->chargeA;
    chargeB = md->chargeB;
    type    = md->typeA;
    typeB   = md->typeB;
    bPert   = md->bPerturbed;
<<<<<<< HEAD
    
    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;
    
    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(cginfo[icg]);
    iwater = GET_CGINFO_SOLOPT(cginfo[icg]);
    
    bFreeEnergy = FALSE;
    if (md->nPerturbed) 
    {
        /* Check if any of the particles involved are perturbed. 
         * If not we can do the cheaper normal put_in_list
         * and use more solvent optimization.
         */
        for(i=0; i<nicg; i++)
=======

    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;

    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(cginfo[icg]);

    iwater = (solvent_opt != esolNO) ? GET_CGINFO_SOLOPT(cginfo[icg]) : esolNO;

    bFreeEnergy = FALSE;
    if (md->nPerturbed)
    {
        /* Check if any of the particles involved are perturbed.
         * If not we can do the cheaper normal put_in_list
         * and use more solvent optimization.
         */
        for (i = 0; i < nicg; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            bFreeEnergy |= bPert[i0+i];
        }
        /* Loop over the j charge groups */
<<<<<<< HEAD
        for(j=0; (j<nj && !bFreeEnergy); j++) 
=======
        for (j = 0; (j < nj && !bFreeEnergy); j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            jcg = jjcg[j];
            jj0 = index[jcg];
            jj1 = index[jcg+1];
<<<<<<< HEAD
            /* Finally loop over the atoms in the j-charge group */	
            for(jj=jj0; jj<jj1; jj++)
=======
            /* Finally loop over the atoms in the j-charge group */
            for (jj = jj0; jj < jj1; jj++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                bFreeEnergy |= bPert[jj];
            }
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Unpack pointers to neighbourlist structs */
    if (fr->nnblists == 1)
    {
        nbl_ind = 0;
    }
    else
    {
<<<<<<< HEAD
        nbl_ind = fr->gid2nblists[GID(igid,jgid,ngid)];
=======
        nbl_ind = fr->gid2nblists[GID(igid, jgid, ngid)];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    if (bLR)
    {
        nlist = fr->nblists[nbl_ind].nlist_lr;
    }
    else
    {
        nlist = fr->nblists[nbl_ind].nlist_sr;
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (iwater != esolNO)
    {
        vdwc = &nlist[eNL_VDWQQ_WATER];
        vdw  = &nlist[eNL_VDW];
        coul = &nlist[eNL_QQ_WATER];
#ifndef DISABLE_WATERWATER_NLIST
        vdwc_ww = &nlist[eNL_VDWQQ_WATERWATER];
        coul_ww = &nlist[eNL_QQ_WATERWATER];
#endif
<<<<<<< HEAD
    } 
    else 
=======
    }
    else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        vdwc = &nlist[eNL_VDWQQ];
        vdw  = &nlist[eNL_VDW];
        coul = &nlist[eNL_QQ];
    }
<<<<<<< HEAD
    
    if (!bFreeEnergy) 
    {
        if (iwater != esolNO) 
        {
            /* Loop over the atoms in the i charge group */    
            i_atom  = i0;
            gid     = GID(igid,jgid,ngid);
            /* Create new i_atom for each energy group */
            if (bDoCoul && bDoVdW)
            {
                new_i_nblist(vdwc,bLR,i_atom,shift,gid);
#ifndef DISABLE_WATERWATER_NLIST
                new_i_nblist(vdwc_ww,bLR,i_atom,shift,gid);
=======

    if (!bFreeEnergy)
    {
        if (iwater != esolNO)
        {
            /* Loop over the atoms in the i charge group */
            i_atom  = i0;
            gid     = GID(igid, jgid, ngid);
            /* Create new i_atom for each energy group */
            if (bDoCoul && bDoVdW)
            {
                new_i_nblist(vdwc, bLR, i_atom, shift, gid);
#ifndef DISABLE_WATERWATER_NLIST
                new_i_nblist(vdwc_ww, bLR, i_atom, shift, gid);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
            }
            if (bDoVdW)
            {
<<<<<<< HEAD
                new_i_nblist(vdw,bLR,i_atom,shift,gid);
            }
            if (bDoCoul) 
            {
                new_i_nblist(coul,bLR,i_atom,shift,gid);
#ifndef DISABLE_WATERWATER_NLIST
                new_i_nblist(coul_ww,bLR,i_atom,shift,gid);
#endif
            }      
	  /* Loop over the j charge groups */
            for(j=0; (j<nj); j++) 
            {
                jcg=jjcg[j];
                
=======
                new_i_nblist(vdw, bLR, i_atom, shift, gid);
            }
            if (bDoCoul)
            {
                new_i_nblist(coul, bLR, i_atom, shift, gid);
#ifndef DISABLE_WATERWATER_NLIST
                new_i_nblist(coul_ww, bLR, i_atom, shift, gid);
#endif
            }
            /* Loop over the j charge groups */
            for (j = 0; (j < nj); j++)
            {
                jcg = jjcg[j];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                if (jcg == icg)
                {
                    continue;
                }
<<<<<<< HEAD
                
                jj0 = index[jcg];
                jwater = GET_CGINFO_SOLOPT(cginfo[jcg]);
                
=======

                jj0    = index[jcg];
                jwater = GET_CGINFO_SOLOPT(cginfo[jcg]);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                if (iwater == esolSPC && jwater == esolSPC)
                {
                    /* Interaction between two SPC molecules */
                    if (!bDoCoul)
                    {
                        /* VdW only - only first atoms in each water interact */
<<<<<<< HEAD
                        add_j_to_nblist(vdw,jj0,bLR);
                    }
                    else 
                    {
#ifdef DISABLE_WATERWATER_NLIST	
                        /* Add entries for the three atoms - only do VdW if we need to */
                        if (!bDoVdW)
                        {
                            add_j_to_nblist(coul,jj0,bLR);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc,jj0,bLR);
                        }
                        add_j_to_nblist(coul,jj0+1,bLR);
                        add_j_to_nblist(coul,jj0+2,bLR);	    
=======
                        add_j_to_nblist(vdw, jj0, bLR);
                    }
                    else
                    {
#ifdef DISABLE_WATERWATER_NLIST
                        /* Add entries for the three atoms - only do VdW if we need to */
                        if (!bDoVdW)
                        {
                            add_j_to_nblist(coul, jj0, bLR);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc, jj0, bLR);
                        }
                        add_j_to_nblist(coul, jj0+1, bLR);
                        add_j_to_nblist(coul, jj0+2, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#else
                        /* One entry for the entire water-water interaction */
                        if (!bDoVdW)
                        {
<<<<<<< HEAD
                            add_j_to_nblist(coul_ww,jj0,bLR);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc_ww,jj0,bLR);
                        }
#endif
                    }  
                } 
                else if (iwater == esolTIP4P && jwater == esolTIP4P) 
=======
                            add_j_to_nblist(coul_ww, jj0, bLR);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc_ww, jj0, bLR);
                        }
#endif
                    }
                }
                else if (iwater == esolTIP4P && jwater == esolTIP4P)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    /* Interaction between two TIP4p molecules */
                    if (!bDoCoul)
                    {
                        /* VdW only - only first atoms in each water interact */
<<<<<<< HEAD
                        add_j_to_nblist(vdw,jj0,bLR);
                    }
                    else 
                    {
#ifdef DISABLE_WATERWATER_NLIST	
                        /* Add entries for the four atoms - only do VdW if we need to */
                        if (bDoVdW)
                        {
                            add_j_to_nblist(vdw,jj0,bLR);
                        }
                        add_j_to_nblist(coul,jj0+1,bLR);
                        add_j_to_nblist(coul,jj0+2,bLR);	    
                        add_j_to_nblist(coul,jj0+3,bLR);	    
=======
                        add_j_to_nblist(vdw, jj0, bLR);
                    }
                    else
                    {
#ifdef DISABLE_WATERWATER_NLIST
                        /* Add entries for the four atoms - only do VdW if we need to */
                        if (bDoVdW)
                        {
                            add_j_to_nblist(vdw, jj0, bLR);
                        }
                        add_j_to_nblist(coul, jj0+1, bLR);
                        add_j_to_nblist(coul, jj0+2, bLR);
                        add_j_to_nblist(coul, jj0+3, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#else
                        /* One entry for the entire water-water interaction */
                        if (!bDoVdW)
                        {
<<<<<<< HEAD
                            add_j_to_nblist(coul_ww,jj0,bLR);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc_ww,jj0,bLR);
                        }
#endif
                    }  					
                }
                else 
=======
                            add_j_to_nblist(coul_ww, jj0, bLR);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc_ww, jj0, bLR);
                        }
#endif
                    }
                }
                else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    /* j charge group is not water, but i is.
                     * Add entries to the water-other_atom lists; the geometry of the water
                     * molecule doesn't matter - that is taken care of in the nonbonded kernel,
                     * so we don't care if it is SPC or TIP4P...
                     */
<<<<<<< HEAD
                    
                    jj1 = index[jcg+1];
                    
                    if (!bDoVdW) 
                    {
                        for(jj=jj0; (jj<jj1); jj++) 
                        {
                            if (charge[jj] != 0)
                            {
                                add_j_to_nblist(coul,jj,bLR);
=======

                    jj1 = index[jcg+1];

                    if (!bDoVdW)
                    {
                        for (jj = jj0; (jj < jj1); jj++)
                        {
                            if (charge[jj] != 0)
                            {
                                add_j_to_nblist(coul, jj, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                            }
                        }
                    }
                    else if (!bDoCoul)
                    {
<<<<<<< HEAD
                        for(jj=jj0; (jj<jj1); jj++)
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                add_j_to_nblist(vdw,jj,bLR);
                            }
                        }
                    }
                    else 
                    {
                        /* _charge_ _groups_ interact with both coulomb and LJ */
                        /* Check which atoms we should add to the lists!       */
                        for(jj=jj0; (jj<jj1); jj++) 
                        {
                            if (bHaveVdW[type[jj]]) 
                            {
                                if (charge[jj] != 0)
                                {
                                    add_j_to_nblist(vdwc,jj,bLR);
                                }
                                else
                                {
                                    add_j_to_nblist(vdw,jj,bLR);
=======
                        for (jj = jj0; (jj < jj1); jj++)
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                add_j_to_nblist(vdw, jj, bLR);
                            }
                        }
                    }
                    else
                    {
                        /* _charge_ _groups_ interact with both coulomb and LJ */
                        /* Check which atoms we should add to the lists!       */
                        for (jj = jj0; (jj < jj1); jj++)
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                if (charge[jj] != 0)
                                {
                                    add_j_to_nblist(vdwc, jj, bLR);
                                }
                                else
                                {
                                    add_j_to_nblist(vdw, jj, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                }
                            }
                            else if (charge[jj] != 0)
                            {
<<<<<<< HEAD
                                add_j_to_nblist(coul,jj,bLR);
=======
                                add_j_to_nblist(coul, jj, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                            }
                        }
                    }
                }
            }
<<<<<<< HEAD
            close_i_nblist(vdw); 
            close_i_nblist(coul); 
            close_i_nblist(vdwc);  
#ifndef DISABLE_WATERWATER_NLIST
            close_i_nblist(coul_ww);
            close_i_nblist(vdwc_ww); 
#endif
        } 
        else
        { 
            /* no solvent as i charge group */
            /* Loop over the atoms in the i charge group */    
            for(i=0; i<nicg; i++) 
            {
                i_atom  = i0+i;
                gid     = GID(igid,jgid,ngid);
                qi      = charge[i_atom];
                
                /* Create new i_atom for each energy group */
                if (bDoVdW && bDoCoul)
                {
                    new_i_nblist(vdwc,bLR,i_atom,shift,gid);
                }
                if (bDoVdW)
                {
                    new_i_nblist(vdw,bLR,i_atom,shift,gid);
                }
                if (bDoCoul)
                {
                    new_i_nblist(coul,bLR,i_atom,shift,gid);
                }
                bDoVdW_i  = (bDoVdW  && bHaveVdW[type[i_atom]]);
                bDoCoul_i = (bDoCoul && qi!=0);
                
                if (bDoVdW_i || bDoCoul_i) 
                {
                    /* Loop over the j charge groups */
                    for(j=0; (j<nj); j++) 
                    {
                        jcg=jjcg[j];
                        
=======
            close_i_nblist(vdw);
            close_i_nblist(coul);
            close_i_nblist(vdwc);
#ifndef DISABLE_WATERWATER_NLIST
            close_i_nblist(coul_ww);
            close_i_nblist(vdwc_ww);
#endif
        }
        else
        {
            /* no solvent as i charge group */
            /* Loop over the atoms in the i charge group */
            for (i = 0; i < nicg; i++)
            {
                i_atom  = i0+i;
                gid     = GID(igid, jgid, ngid);
                qi      = charge[i_atom];

                /* Create new i_atom for each energy group */
                if (bDoVdW && bDoCoul)
                {
                    new_i_nblist(vdwc, bLR, i_atom, shift, gid);
                }
                if (bDoVdW)
                {
                    new_i_nblist(vdw, bLR, i_atom, shift, gid);
                }
                if (bDoCoul)
                {
                    new_i_nblist(coul, bLR, i_atom, shift, gid);
                }
                bDoVdW_i  = (bDoVdW  && bHaveVdW[type[i_atom]]);
                bDoCoul_i = (bDoCoul && qi != 0);

                if (bDoVdW_i || bDoCoul_i)
                {
                    /* Loop over the j charge groups */
                    for (j = 0; (j < nj); j++)
                    {
                        jcg = jjcg[j];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                        /* Check for large charge groups */
                        if (jcg == icg)
                        {
                            jj0 = i0 + i + 1;
                        }
                        else
                        {
                            jj0 = index[jcg];
                        }
<<<<<<< HEAD
                        
                        jj1=index[jcg+1];
                        /* Finally loop over the atoms in the j-charge group */	
                        for(jj=jj0; jj<jj1; jj++) 
                        {
                            bNotEx = NOTEXCL(bExcl,i,jj);
                            
                            if (bNotEx) 
                            {
                                if (!bDoVdW_i) 
                                { 
                                    if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul,jj,bLR);
                                    }
                                }
                                else if (!bDoCoul_i) 
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        add_j_to_nblist(vdw,jj,bLR);
                                    }
                                }
                                else 
                                {
                                    if (bHaveVdW[type[jj]]) 
                                    {
                                        if (charge[jj] != 0)
                                        {
                                            add_j_to_nblist(vdwc,jj,bLR);
                                        }
                                        else
                                        {
                                            add_j_to_nblist(vdw,jj,bLR);
                                        }
                                    } 
                                    else if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul,jj,bLR);
=======

                        jj1 = index[jcg+1];
                        /* Finally loop over the atoms in the j-charge group */
                        for (jj = jj0; jj < jj1; jj++)
                        {
                            bNotEx = NOTEXCL(bExcl, i, jj);

                            if (bNotEx)
                            {
                                if (!bDoVdW_i)
                                {
                                    if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj, bLR);
                                    }
                                }
                                else if (!bDoCoul_i)
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        add_j_to_nblist(vdw, jj, bLR);
                                    }
                                }
                                else
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        if (charge[jj] != 0)
                                        {
                                            add_j_to_nblist(vdwc, jj, bLR);
                                        }
                                        else
                                        {
                                            add_j_to_nblist(vdw, jj, bLR);
                                        }
                                    }
                                    else if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                    }
                                }
                            }
                        }
                    }
                }
                close_i_nblist(vdw);
                close_i_nblist(coul);
                close_i_nblist(vdwc);
            }
        }
    }
    else
    {
        /* we are doing free energy */
        vdwc_free = &nlist[eNL_VDWQQ_FREE];
        vdw_free  = &nlist[eNL_VDW_FREE];
        coul_free = &nlist[eNL_QQ_FREE];
<<<<<<< HEAD
        /* Loop over the atoms in the i charge group */    
        for(i=0; i<nicg; i++) 
        {
            i_atom  = i0+i;
            gid     = GID(igid,jgid,ngid);
            qi      = charge[i_atom];
            qiB     = chargeB[i_atom];
            
            /* Create new i_atom for each energy group */
            if (bDoVdW && bDoCoul) 
                new_i_nblist(vdwc,bLR,i_atom,shift,gid);
            if (bDoVdW)   
                new_i_nblist(vdw,bLR,i_atom,shift,gid);
            if (bDoCoul) 
                new_i_nblist(coul,bLR,i_atom,shift,gid);
            
            new_i_nblist(vdw_free,bLR,i_atom,shift,gid);
            new_i_nblist(coul_free,bLR,i_atom,shift,gid);
            new_i_nblist(vdwc_free,bLR,i_atom,shift,gid);
            
            bDoVdW_i  = (bDoVdW  &&
                         (bHaveVdW[type[i_atom]] || bHaveVdW[typeB[i_atom]]));
            bDoCoul_i = (bDoCoul && (qi!=0 || qiB!=0));
=======
        /* Loop over the atoms in the i charge group */
        for (i = 0; i < nicg; i++)
        {
            i_atom  = i0+i;
            gid     = GID(igid, jgid, ngid);
            qi      = charge[i_atom];
            qiB     = chargeB[i_atom];

            /* Create new i_atom for each energy group */
            if (bDoVdW && bDoCoul)
            {
                new_i_nblist(vdwc, bLR, i_atom, shift, gid);
            }
            if (bDoVdW)
            {
                new_i_nblist(vdw, bLR, i_atom, shift, gid);
            }
            if (bDoCoul)
            {
                new_i_nblist(coul, bLR, i_atom, shift, gid);
            }

            new_i_nblist(vdw_free, bLR, i_atom, shift, gid);
            new_i_nblist(coul_free, bLR, i_atom, shift, gid);
            new_i_nblist(vdwc_free, bLR, i_atom, shift, gid);

            bDoVdW_i  = (bDoVdW  &&
                         (bHaveVdW[type[i_atom]] || bHaveVdW[typeB[i_atom]]));
            bDoCoul_i = (bDoCoul && (qi != 0 || qiB != 0));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* For TIP4P the first atom does not have a charge,
             * but the last three do. So we should still put an atom
             * without LJ but with charge in the water-atom neighborlist
             * for a TIP4p i charge group.
             * For SPC type water the first atom has LJ and charge,
             * so there is no such problem.
             */
            if (iwater == esolNO)
            {
                bDoCoul_i_sol = bDoCoul_i;
            }
            else
            {
                bDoCoul_i_sol = bDoCoul;
            }
<<<<<<< HEAD
            
            if (bDoVdW_i || bDoCoul_i_sol) 
            {
                /* Loop over the j charge groups */
                for(j=0; (j<nj); j++)
                {
                    jcg=jjcg[j];
                    
=======

            if (bDoVdW_i || bDoCoul_i_sol)
            {
                /* Loop over the j charge groups */
                for (j = 0; (j < nj); j++)
                {
                    jcg = jjcg[j];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    /* Check for large charge groups */
                    if (jcg == icg)
                    {
                        jj0 = i0 + i + 1;
                    }
                    else
                    {
                        jj0 = index[jcg];
                    }
<<<<<<< HEAD
                    
                    jj1=index[jcg+1];
                    /* Finally loop over the atoms in the j-charge group */	
                    bFree = bPert[i_atom];
                    for(jj=jj0; (jj<jj1); jj++) 
=======

                    jj1 = index[jcg+1];
                    /* Finally loop over the atoms in the j-charge group */
                    bFree = bPert[i_atom];
                    for (jj = jj0; (jj < jj1); jj++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        bFreeJ = bFree || bPert[jj];
                        /* Complicated if, because the water H's should also
                         * see perturbed j-particles
                         */
<<<<<<< HEAD
                        if (iwater==esolNO || i==0 || bFreeJ) 
                        {
                            bNotEx = NOTEXCL(bExcl,i,jj);
                            
                            if (bNotEx) 
                            {
                                if (bFreeJ)
                                {
                                    if (!bDoVdW_i) 
                                    {
                                        if (charge[jj]!=0 || chargeB[jj]!=0)
                                        {
                                            add_j_to_nblist(coul_free,jj,bLR);
                                        }
                                    }
                                    else if (!bDoCoul_i) 
                                    {
                                        if (bHaveVdW[type[jj]] || bHaveVdW[typeB[jj]])
                                        {
                                            add_j_to_nblist(vdw_free,jj,bLR);
                                        }
                                    }
                                    else 
                                    {
                                        if (bHaveVdW[type[jj]] || bHaveVdW[typeB[jj]]) 
                                        {
                                            if (charge[jj]!=0 || chargeB[jj]!=0)
                                            {
                                                add_j_to_nblist(vdwc_free,jj,bLR);
                                            }
                                            else
                                            {
                                                add_j_to_nblist(vdw_free,jj,bLR);
                                            }
                                        }
                                        else if (charge[jj]!=0 || chargeB[jj]!=0)
                                            add_j_to_nblist(coul_free,jj,bLR);
                                    }
                                }
                                else if (!bDoVdW_i) 
                                { 
                                    /* This is done whether or not bWater is set */
                                    if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul,jj,bLR);
                                    }
                                }
                                else if (!bDoCoul_i_sol) 
                                { 
                                    if (bHaveVdW[type[jj]])
                                    {
                                        add_j_to_nblist(vdw,jj,bLR);
                                    }
                                }
                                else 
                                {
                                    if (bHaveVdW[type[jj]]) 
                                    {
                                        if (charge[jj] != 0)
                                        {
                                            add_j_to_nblist(vdwc,jj,bLR);
                                        }
                                        else
                                        {
                                            add_j_to_nblist(vdw,jj,bLR);
                                        }
                                    } 
                                    else if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul,jj,bLR);
=======
                        if (iwater == esolNO || i == 0 || bFreeJ)
                        {
                            bNotEx = NOTEXCL(bExcl, i, jj);

                            if (bNotEx)
                            {
                                if (bFreeJ)
                                {
                                    if (!bDoVdW_i)
                                    {
                                        if (charge[jj] != 0 || chargeB[jj] != 0)
                                        {
                                            add_j_to_nblist(coul_free, jj, bLR);
                                        }
                                    }
                                    else if (!bDoCoul_i)
                                    {
                                        if (bHaveVdW[type[jj]] || bHaveVdW[typeB[jj]])
                                        {
                                            add_j_to_nblist(vdw_free, jj, bLR);
                                        }
                                    }
                                    else
                                    {
                                        if (bHaveVdW[type[jj]] || bHaveVdW[typeB[jj]])
                                        {
                                            if (charge[jj] != 0 || chargeB[jj] != 0)
                                            {
                                                add_j_to_nblist(vdwc_free, jj, bLR);
                                            }
                                            else
                                            {
                                                add_j_to_nblist(vdw_free, jj, bLR);
                                            }
                                        }
                                        else if (charge[jj] != 0 || chargeB[jj] != 0)
                                        {
                                            add_j_to_nblist(coul_free, jj, bLR);
                                        }
                                    }
                                }
                                else if (!bDoVdW_i)
                                {
                                    /* This is done whether or not bWater is set */
                                    if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj, bLR);
                                    }
                                }
                                else if (!bDoCoul_i_sol)
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        add_j_to_nblist(vdw, jj, bLR);
                                    }
                                }
                                else
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        if (charge[jj] != 0)
                                        {
                                            add_j_to_nblist(vdwc, jj, bLR);
                                        }
                                        else
                                        {
                                            add_j_to_nblist(vdw, jj, bLR);
                                        }
                                    }
                                    else if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj, bLR);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                    }
                                }
                            }
                        }
                    }
                }
            }
            close_i_nblist(vdw);
            close_i_nblist(coul);
            close_i_nblist(vdwc);
            close_i_nblist(vdw_free);
            close_i_nblist(coul_free);
            close_i_nblist(vdwc_free);
        }
    }
}

<<<<<<< HEAD
static void 
put_in_list_qmmm(gmx_bool              bHaveVdW[],
                 int               ngid,
                 t_mdatoms *       md,
                 int               icg,
                 int               jgid,
                 int               nj,
                 atom_id           jjcg[],
                 atom_id           index[],
                 t_excl            bExcl[],
                 int               shift,
                 t_forcerec *      fr,
                 gmx_bool              bLR,
                 gmx_bool              bDoVdW,
                 gmx_bool              bDoCoul)
{
    t_nblist *   coul;
    int 	  i,j,jcg,igid,gid;
    atom_id   jj,jj0,jj1,i_atom;
    int       i0,nicg;
    gmx_bool      bNotEx;
    
    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;
    
    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(fr->cginfo[icg]);
    
    coul = &fr->QMMMlist;
    
    /* Loop over atoms in the ith charge group */
    for (i=0;i<nicg;i++)
    {
        i_atom = i0+i;
        gid    = GID(igid,jgid,ngid);
        /* Create new i_atom for each energy group */
        new_i_nblist(coul,bLR,i_atom,shift,gid);
        
        /* Loop over the j charge groups */
        for (j=0;j<nj;j++)
        {
            jcg=jjcg[j];
            
            /* Charge groups cannot have QM and MM atoms simultaneously */
            if (jcg!=icg)
=======
static void
put_in_list_adress(gmx_bool              bHaveVdW[],
                   int                   ngid,
                   t_mdatoms     *       md,
                   int                   icg,
                   int                   jgid,
                   int                   nj,
                   atom_id               jjcg[],
                   atom_id               index[],
                   t_excl                bExcl[],
                   int                   shift,
                   t_forcerec     *      fr,
                   gmx_bool              bLR,
                   gmx_bool              bDoVdW,
                   gmx_bool              bDoCoul,
                   int                   solvent_opt)
{
    /* The a[] index has been removed,
     * to put it back in i_atom should be a[i0] and jj should be a[jj].
     */
    t_nblist  *   vdwc;
    t_nblist  *   vdw;
    t_nblist  *   coul;
    t_nblist  *   vdwc_adress  = NULL;
    t_nblist  *   vdw_adress   = NULL;
    t_nblist  *   coul_adress  = NULL;
    t_nblist  *   vdwc_ww      = NULL;
    t_nblist  *   coul_ww      = NULL;

    int           i, j, jcg, igid, gid, nbl_ind, nbl_ind_adress;
    atom_id       jj, jj0, jj1, i_atom;
    int           i0, nicg, len;

    int          *cginfo;
    int          *type, *typeB;
    real         *charge, *chargeB;
    real         *wf;
    real          qi, qiB, qq, rlj;
    gmx_bool      bFreeEnergy, bFree, bFreeJ, bNotEx, *bPert;
    gmx_bool      bDoVdW_i, bDoCoul_i, bDoCoul_i_sol;
    gmx_bool      b_hybrid;
    gmx_bool      j_all_atom;
    int           iwater, jwater;
    t_nblist     *nlist, *nlist_adress;
    gmx_bool      bEnergyGroupCG;

    /* Copy some pointers */
    cginfo  = fr->cginfo;
    charge  = md->chargeA;
    chargeB = md->chargeB;
    type    = md->typeA;
    typeB   = md->typeB;
    bPert   = md->bPerturbed;
    wf      = md->wf;

    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;

    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(cginfo[icg]);

    iwater = (solvent_opt != esolNO) ? GET_CGINFO_SOLOPT(cginfo[icg]) : esolNO;

    if (md->nPerturbed)
    {
        gmx_fatal(FARGS, "AdResS does not support free energy pertubation\n");
    }

    /* Unpack pointers to neighbourlist structs */
    if (fr->nnblists == 2)
    {
        nbl_ind        = 0;
        nbl_ind_adress = 1;
    }
    else
    {
        nbl_ind        = fr->gid2nblists[GID(igid, jgid, ngid)];
        nbl_ind_adress = nbl_ind+fr->nnblists/2;
    }
    if (bLR)
    {
        nlist        = fr->nblists[nbl_ind].nlist_lr;
        nlist_adress = fr->nblists[nbl_ind_adress].nlist_lr;
    }
    else
    {
        nlist        = fr->nblists[nbl_ind].nlist_sr;
        nlist_adress = fr->nblists[nbl_ind_adress].nlist_sr;
    }


    vdwc = &nlist[eNL_VDWQQ];
    vdw  = &nlist[eNL_VDW];
    coul = &nlist[eNL_QQ];

    vdwc_adress = &nlist_adress[eNL_VDWQQ];
    vdw_adress  = &nlist_adress[eNL_VDW];
    coul_adress = &nlist_adress[eNL_QQ];

    /* We do not support solvent optimization with AdResS for now.
       For this we would need hybrid solvent-other kernels */

    /* no solvent as i charge group */
    /* Loop over the atoms in the i charge group */
    for (i = 0; i < nicg; i++)
    {
        i_atom  = i0+i;
        gid     = GID(igid, jgid, ngid);
        qi      = charge[i_atom];

        /* Create new i_atom for each energy group */
        if (bDoVdW && bDoCoul)
        {
            new_i_nblist(vdwc, bLR, i_atom, shift, gid);
            new_i_nblist(vdwc_adress, bLR, i_atom, shift, gid);

        }
        if (bDoVdW)
        {
            new_i_nblist(vdw, bLR, i_atom, shift, gid);
            new_i_nblist(vdw_adress, bLR, i_atom, shift, gid);

        }
        if (bDoCoul)
        {
            new_i_nblist(coul, bLR, i_atom, shift, gid);
            new_i_nblist(coul_adress, bLR, i_atom, shift, gid);
        }
        bDoVdW_i  = (bDoVdW  && bHaveVdW[type[i_atom]]);
        bDoCoul_i = (bDoCoul && qi != 0);

        /* Here we find out whether the energy groups interaction belong to a
         * coarse-grained (vsite) or atomistic interaction. Note that, beacuse
         * interactions between coarse-grained and other (atomistic) energygroups
         * are excluded automatically by grompp, it is sufficient to check for
         * the group id of atom i (igid) */
        bEnergyGroupCG = !egp_explicit(fr, igid);

        if (bDoVdW_i || bDoCoul_i)
        {
            /* Loop over the j charge groups */
            for (j = 0; (j < nj); j++)
            {
                jcg = jjcg[j];

                /* Check for large charge groups */
                if (jcg == icg)
                {
                    jj0 = i0 + i + 1;
                }
                else
                {
                    jj0 = index[jcg];
                }

                jj1 = index[jcg+1];
                /* Finally loop over the atoms in the j-charge group */
                for (jj = jj0; jj < jj1; jj++)
                {
                    bNotEx = NOTEXCL(bExcl, i, jj);

                    /* Now we have to exclude interactions which will be zero
                     * anyway due to the AdResS weights (in previous implementations
                     * this was done in the force kernel). This is necessary as
                     * pure interactions (those with b_hybrid=false, i.e. w_i*w_j==1 or 0)
                     * are put into neighbour lists which will be passed to the
                     * standard (optimized) kernels for speed. The interactions with
                     * b_hybrid=true are placed into the _adress neighbour lists and
                     * processed by the generic AdResS kernel.
                     */
                    if ( (bEnergyGroupCG &&
                          wf[i_atom] >= 1-GMX_REAL_EPS && wf[jj] >= 1-GMX_REAL_EPS ) ||
                         ( !bEnergyGroupCG && wf[jj] <= GMX_REAL_EPS ) )
                    {
                        continue;
                    }

                    b_hybrid = !((wf[i_atom] >= 1-GMX_REAL_EPS && wf[jj] >= 1-GMX_REAL_EPS) ||
                                 (wf[i_atom] <= GMX_REAL_EPS && wf[jj] <= GMX_REAL_EPS));

                    if (bNotEx)
                    {
                        if (!bDoVdW_i)
                        {
                            if (charge[jj] != 0)
                            {
                                if (!b_hybrid)
                                {
                                    add_j_to_nblist(coul, jj, bLR);
                                }
                                else
                                {
                                    add_j_to_nblist(coul_adress, jj, bLR);
                                }
                            }
                        }
                        else if (!bDoCoul_i)
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                if (!b_hybrid)
                                {
                                    add_j_to_nblist(vdw, jj, bLR);
                                }
                                else
                                {
                                    add_j_to_nblist(vdw_adress, jj, bLR);
                                }
                            }
                        }
                        else
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                if (charge[jj] != 0)
                                {
                                    if (!b_hybrid)
                                    {
                                        add_j_to_nblist(vdwc, jj, bLR);
                                    }
                                    else
                                    {
                                        add_j_to_nblist(vdwc_adress, jj, bLR);
                                    }
                                }
                                else
                                {
                                    if (!b_hybrid)
                                    {
                                        add_j_to_nblist(vdw, jj, bLR);
                                    }
                                    else
                                    {
                                        add_j_to_nblist(vdw_adress, jj, bLR);
                                    }

                                }
                            }
                            else if (charge[jj] != 0)
                            {
                                if (!b_hybrid)
                                {
                                    add_j_to_nblist(coul, jj, bLR);
                                }
                                else
                                {
                                    add_j_to_nblist(coul_adress, jj, bLR);
                                }

                            }
                        }
                    }
                }
            }

            close_i_nblist(vdw);
            close_i_nblist(coul);
            close_i_nblist(vdwc);
            close_i_nblist(vdw_adress);
            close_i_nblist(coul_adress);
            close_i_nblist(vdwc_adress);
        }
    }
}

static void
put_in_list_qmmm(gmx_bool              bHaveVdW[],
                 int                   ngid,
                 t_mdatoms     *       md,
                 int                   icg,
                 int                   jgid,
                 int                   nj,
                 atom_id               jjcg[],
                 atom_id               index[],
                 t_excl                bExcl[],
                 int                   shift,
                 t_forcerec     *      fr,
                 gmx_bool              bLR,
                 gmx_bool              bDoVdW,
                 gmx_bool              bDoCoul,
                 int                   solvent_opt)
{
    t_nblist  *   coul;
    int           i, j, jcg, igid, gid;
    atom_id       jj, jj0, jj1, i_atom;
    int           i0, nicg;
    gmx_bool      bNotEx;

    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;

    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(fr->cginfo[icg]);

    coul = &fr->QMMMlist;

    /* Loop over atoms in the ith charge group */
    for (i = 0; i < nicg; i++)
    {
        i_atom = i0+i;
        gid    = GID(igid, jgid, ngid);
        /* Create new i_atom for each energy group */
        new_i_nblist(coul, bLR, i_atom, shift, gid);

        /* Loop over the j charge groups */
        for (j = 0; j < nj; j++)
        {
            jcg = jjcg[j];

            /* Charge groups cannot have QM and MM atoms simultaneously */
            if (jcg != icg)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                jj0 = index[jcg];
                jj1 = index[jcg+1];
                /* Finally loop over the atoms in the j-charge group */
<<<<<<< HEAD
                for(jj=jj0; jj<jj1; jj++)
                {
                    bNotEx = NOTEXCL(bExcl,i,jj);
                    if(bNotEx)
                        add_j_to_nblist(coul,jj,bLR);
=======
                for (jj = jj0; jj < jj1; jj++)
                {
                    bNotEx = NOTEXCL(bExcl, i, jj);
                    if (bNotEx)
                    {
                        add_j_to_nblist(coul, jj, bLR);
                    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
            }
        }
        close_i_nblist(coul);
    }
}

<<<<<<< HEAD
static void 
put_in_list_cg(gmx_bool              bHaveVdW[],
               int               ngid,
               t_mdatoms *       md,
               int               icg,
               int               jgid,
               int               nj,
               atom_id           jjcg[],
               atom_id           index[],
               t_excl            bExcl[],
               int               shift,
               t_forcerec *      fr,
               gmx_bool              bLR,
               gmx_bool              bDoVdW,
               gmx_bool              bDoCoul)
{
    int          cginfo;
    int          igid,gid,nbl_ind;
    t_nblist *   vdwc;
    int          j,jcg;
=======
static void
put_in_list_cg(gmx_bool              bHaveVdW[],
               int                   ngid,
               t_mdatoms     *       md,
               int                   icg,
               int                   jgid,
               int                   nj,
               atom_id               jjcg[],
               atom_id               index[],
               t_excl                bExcl[],
               int                   shift,
               t_forcerec     *      fr,
               gmx_bool              bLR,
               gmx_bool              bDoVdW,
               gmx_bool              bDoCoul,
               int                   solvent_opt)
{
    int          cginfo;
    int          igid, gid, nbl_ind;
    t_nblist *   vdwc;
    int          j, jcg;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    cginfo = fr->cginfo[icg];

    igid = GET_CGINFO_GID(cginfo);
<<<<<<< HEAD
    gid  = GID(igid,jgid,ngid);
=======
    gid  = GID(igid, jgid, ngid);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Unpack pointers to neighbourlist structs */
    if (fr->nnblists == 1)
    {
        nbl_ind = 0;
    }
    else
    {
        nbl_ind = fr->gid2nblists[gid];
    }
    if (bLR)
    {
        vdwc = &fr->nblists[nbl_ind].nlist_lr[eNL_VDWQQ];
    }
    else
    {
        vdwc = &fr->nblists[nbl_ind].nlist_sr[eNL_VDWQQ];
    }

    /* Make a new neighbor list for charge group icg.
     * Currently simply one neighbor list is made with LJ and Coulomb.
     * If required, zero interactions could be removed here
     * or in the force loop.
     */
<<<<<<< HEAD
    new_i_nblist(vdwc,bLR,index[icg],shift,gid);
    vdwc->iinr_end[vdwc->nri] = index[icg+1];

    for(j=0; (j<nj); j++) 
=======
    new_i_nblist(vdwc, bLR, index[icg], shift, gid);
    vdwc->iinr_end[vdwc->nri] = index[icg+1];

    for (j = 0; (j < nj); j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        jcg = jjcg[j];
        /* Skip the icg-icg pairs if all self interactions are excluded */
        if (!(jcg == icg && GET_CGINFO_EXCL_INTRA(cginfo)))
        {
            /* Here we add the j charge group jcg to the list,
             * exclusions are also added to the list.
             */
<<<<<<< HEAD
            add_j_to_nblist_cg(vdwc,index[jcg],index[jcg+1],bExcl,icg==jcg,bLR);
        }
    }

    close_i_nblist(vdwc);  
}

static void setexcl(atom_id start,atom_id end,t_blocka *excl,gmx_bool b,
                    t_excl bexcl[])
{
    atom_id i,k;
    
    if (b)
    {
        for(i=start; i<end; i++)
        {
            for(k=excl->index[i]; k<excl->index[i+1]; k++)
            {
                SETEXCL(bexcl,i-start,excl->a[k]);
=======
            add_j_to_nblist_cg(vdwc, index[jcg], index[jcg+1], bExcl, icg == jcg, bLR);
        }
    }

    close_i_nblist(vdwc);
}

static void setexcl(atom_id start, atom_id end, t_blocka *excl, gmx_bool b,
                    t_excl bexcl[])
{
    atom_id i, k;

    if (b)
    {
        for (i = start; i < end; i++)
        {
            for (k = excl->index[i]; k < excl->index[i+1]; k++)
            {
                SETEXCL(bexcl, i-start, excl->a[k]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
    else
    {
<<<<<<< HEAD
        for(i=start; i<end; i++)
        {
            for(k=excl->index[i]; k<excl->index[i+1]; k++)
            {
                RMEXCL(bexcl,i-start,excl->a[k]);
=======
        for (i = start; i < end; i++)
        {
            for (k = excl->index[i]; k < excl->index[i+1]; k++)
            {
                RMEXCL(bexcl, i-start, excl->a[k]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
}

<<<<<<< HEAD
int calc_naaj(int icg,int cgtot)
{
    int naaj;
    
=======
int calc_naaj(int icg, int cgtot)
{
    int naaj;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if ((cgtot % 2) == 1)
    {
        /* Odd number of charge groups, easy */
        naaj = 1 + (cgtot/2);
    }
    else if ((cgtot % 4) == 0)
    {
<<<<<<< HEAD
    /* Multiple of four is hard */
=======
        /* Multiple of four is hard */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (icg < cgtot/2)
        {
            if ((icg % 2) == 0)
            {
<<<<<<< HEAD
                naaj=1+(cgtot/2);
            }
            else
            {
                naaj=cgtot/2;
=======
                naaj = 1+(cgtot/2);
            }
            else
            {
                naaj = cgtot/2;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
        else
        {
            if ((icg % 2) == 1)
            {
<<<<<<< HEAD
                naaj=1+(cgtot/2);
            }
            else
            {
                naaj=cgtot/2;
=======
                naaj = 1+(cgtot/2);
            }
            else
            {
                naaj = cgtot/2;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
    else
    {
        /* cgtot/2 = odd */
        if ((icg % 2) == 0)
        {
<<<<<<< HEAD
            naaj=1+(cgtot/2);
        }
        else
        {
            naaj=cgtot/2;
        }
    }
#ifdef DEBUG
    fprintf(log,"naaj=%d\n",naaj);
=======
            naaj = 1+(cgtot/2);
        }
        else
        {
            naaj = cgtot/2;
        }
    }
#ifdef DEBUG
    fprintf(log, "naaj=%d\n", naaj);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

    return naaj;
}

/************************************************
 *
 *  S I M P L E      C O R E     S T U F F
 *
 ************************************************/

<<<<<<< HEAD
static real calc_image_tric(rvec xi,rvec xj,matrix box,
                            rvec b_inv,int *shift)
=======
static real calc_image_tric(rvec xi, rvec xj, matrix box,
                            rvec b_inv, int *shift)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    /* This code assumes that the cut-off is smaller than
     * a half times the smallest diagonal element of the box.
     */
<<<<<<< HEAD
    const real h25=2.5;
    real dx,dy,dz;
    real r2;
    int  tx,ty,tz;
    
=======
    const real h25 = 2.5;
    real       dx, dy, dz;
    real       r2;
    int        tx, ty, tz;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Compute diff vector */
    dz = xj[ZZ] - xi[ZZ];
    dy = xj[YY] - xi[YY];
    dx = xj[XX] - xi[XX];
<<<<<<< HEAD
    
  /* Perform NINT operation, using trunc operation, therefore
   * we first add 2.5 then subtract 2 again
   */
    tz = dz*b_inv[ZZ] + h25;
=======

    /* Perform NINT operation, using trunc operation, therefore
     * we first add 2.5 then subtract 2 again
     */
    tz  = dz*b_inv[ZZ] + h25;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tz -= 2;
    dz -= tz*box[ZZ][ZZ];
    dy -= tz*box[ZZ][YY];
    dx -= tz*box[ZZ][XX];

<<<<<<< HEAD
    ty = dy*b_inv[YY] + h25;
    ty -= 2;
    dy -= ty*box[YY][YY];
    dx -= ty*box[YY][XX];
    
    tx = dx*b_inv[XX]+h25;
    tx -= 2;
    dx -= tx*box[XX][XX];
  
    /* Distance squared */
    r2 = (dx*dx) + (dy*dy) + (dz*dz);

    *shift = XYZ2IS(tx,ty,tz);
=======
    ty  = dy*b_inv[YY] + h25;
    ty -= 2;
    dy -= ty*box[YY][YY];
    dx -= ty*box[YY][XX];

    tx  = dx*b_inv[XX]+h25;
    tx -= 2;
    dx -= tx*box[XX][XX];

    /* Distance squared */
    r2 = (dx*dx) + (dy*dy) + (dz*dz);

    *shift = XYZ2IS(tx, ty, tz);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return r2;
}

<<<<<<< HEAD
static real calc_image_rect(rvec xi,rvec xj,rvec box_size,
                            rvec b_inv,int *shift)
{
    const real h15=1.5;
    real ddx,ddy,ddz;
    real dx,dy,dz;
    real r2;
    int  tx,ty,tz;
    
=======
static real calc_image_rect(rvec xi, rvec xj, rvec box_size,
                            rvec b_inv, int *shift)
{
    const real h15 = 1.5;
    real       ddx, ddy, ddz;
    real       dx, dy, dz;
    real       r2;
    int        tx, ty, tz;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Compute diff vector */
    dx = xj[XX] - xi[XX];
    dy = xj[YY] - xi[YY];
    dz = xj[ZZ] - xi[ZZ];
<<<<<<< HEAD
  
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Perform NINT operation, using trunc operation, therefore
     * we first add 1.5 then subtract 1 again
     */
    tx = dx*b_inv[XX] + h15;
    ty = dy*b_inv[YY] + h15;
    tz = dz*b_inv[ZZ] + h15;
    tx--;
    ty--;
    tz--;
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Correct diff vector for translation */
    ddx = tx*box_size[XX] - dx;
    ddy = ty*box_size[YY] - dy;
    ddz = tz*box_size[ZZ] - dz;
<<<<<<< HEAD
    
    /* Distance squared */
    r2 = (ddx*ddx) + (ddy*ddy) + (ddz*ddz);
    
    *shift = XYZ2IS(tx,ty,tz);
    
    return r2;
}

static void add_simple(t_ns_buf *nsbuf,int nrj,atom_id cg_j,
                       gmx_bool bHaveVdW[],int ngid,t_mdatoms *md,
                       int icg,int jgid,t_block *cgs,t_excl bexcl[],
                       int shift,t_forcerec *fr,put_in_list_t *put_in_list)
{
    if (nsbuf->nj + nrj > MAX_CG)
    {
        put_in_list(bHaveVdW,ngid,md,icg,jgid,nsbuf->ncg,nsbuf->jcg,
                    cgs->index,bexcl,shift,fr,FALSE,TRUE,TRUE);
=======

    /* Distance squared */
    r2 = (ddx*ddx) + (ddy*ddy) + (ddz*ddz);

    *shift = XYZ2IS(tx, ty, tz);

    return r2;
}

static void add_simple(t_ns_buf *nsbuf, int nrj, atom_id cg_j,
                       gmx_bool bHaveVdW[], int ngid, t_mdatoms *md,
                       int icg, int jgid, t_block *cgs, t_excl bexcl[],
                       int shift, t_forcerec *fr, put_in_list_t *put_in_list)
{
    if (nsbuf->nj + nrj > MAX_CG)
    {
        put_in_list(bHaveVdW, ngid, md, icg, jgid, nsbuf->ncg, nsbuf->jcg,
                    cgs->index, bexcl, shift, fr, FALSE, TRUE, TRUE, fr->solvent_opt);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Reset buffer contents */
        nsbuf->ncg = nsbuf->nj = 0;
    }
    nsbuf->jcg[nsbuf->ncg++] = cg_j;
<<<<<<< HEAD
    nsbuf->nj += nrj;
}

static void ns_inner_tric(rvec x[],int icg,int *i_egp_flags,
                          int njcg,atom_id jcg[],
                          matrix box,rvec b_inv,real rcut2,
                          t_block *cgs,t_ns_buf **ns_buf,
                          gmx_bool bHaveVdW[],int ngid,t_mdatoms *md,
                          t_excl bexcl[],t_forcerec *fr,
                          put_in_list_t *put_in_list)
{
    int      shift;
    int      j,nrj,jgid;
    int      *cginfo=fr->cginfo;
    atom_id  cg_j,*cgindex;
    t_ns_buf *nsbuf;
    
    cgindex = cgs->index;
    shift   = CENTRAL;
    for(j=0; (j<njcg); j++)
    {
        cg_j   = jcg[j];
        nrj    = cgindex[cg_j+1]-cgindex[cg_j];
        if (calc_image_tric(x[icg],x[cg_j],box,b_inv,&shift) < rcut2)
=======
    nsbuf->nj               += nrj;
}

static void ns_inner_tric(rvec x[], int icg, int *i_egp_flags,
                          int njcg, atom_id jcg[],
                          matrix box, rvec b_inv, real rcut2,
                          t_block *cgs, t_ns_buf **ns_buf,
                          gmx_bool bHaveVdW[], int ngid, t_mdatoms *md,
                          t_excl bexcl[], t_forcerec *fr,
                          put_in_list_t *put_in_list)
{
    int       shift;
    int       j, nrj, jgid;
    int      *cginfo = fr->cginfo;
    atom_id   cg_j, *cgindex;
    t_ns_buf *nsbuf;

    cgindex = cgs->index;
    shift   = CENTRAL;
    for (j = 0; (j < njcg); j++)
    {
        cg_j   = jcg[j];
        nrj    = cgindex[cg_j+1]-cgindex[cg_j];
        if (calc_image_tric(x[icg], x[cg_j], box, b_inv, &shift) < rcut2)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            jgid  = GET_CGINFO_GID(cginfo[cg_j]);
            if (!(i_egp_flags[jgid] & EGP_EXCL))
            {
<<<<<<< HEAD
                add_simple(&ns_buf[jgid][shift],nrj,cg_j,
                           bHaveVdW,ngid,md,icg,jgid,cgs,bexcl,shift,fr,
=======
                add_simple(&ns_buf[jgid][shift], nrj, cg_j,
                           bHaveVdW, ngid, md, icg, jgid, cgs, bexcl, shift, fr,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                           put_in_list);
            }
        }
    }
}

<<<<<<< HEAD
static void ns_inner_rect(rvec x[],int icg,int *i_egp_flags,
                          int njcg,atom_id jcg[],
                          gmx_bool bBox,rvec box_size,rvec b_inv,real rcut2,
                          t_block *cgs,t_ns_buf **ns_buf,
                          gmx_bool bHaveVdW[],int ngid,t_mdatoms *md,
                          t_excl bexcl[],t_forcerec *fr,
                          put_in_list_t *put_in_list)
{
    int      shift;
    int      j,nrj,jgid;
    int      *cginfo=fr->cginfo;
    atom_id  cg_j,*cgindex;
=======
static void ns_inner_rect(rvec x[], int icg, int *i_egp_flags,
                          int njcg, atom_id jcg[],
                          gmx_bool bBox, rvec box_size, rvec b_inv, real rcut2,
                          t_block *cgs, t_ns_buf **ns_buf,
                          gmx_bool bHaveVdW[], int ngid, t_mdatoms *md,
                          t_excl bexcl[], t_forcerec *fr,
                          put_in_list_t *put_in_list)
{
    int       shift;
    int       j, nrj, jgid;
    int      *cginfo = fr->cginfo;
    atom_id   cg_j, *cgindex;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_ns_buf *nsbuf;

    cgindex = cgs->index;
    if (bBox)
    {
        shift = CENTRAL;
<<<<<<< HEAD
        for(j=0; (j<njcg); j++)
        {
            cg_j   = jcg[j];
            nrj    = cgindex[cg_j+1]-cgindex[cg_j];
            if (calc_image_rect(x[icg],x[cg_j],box_size,b_inv,&shift) < rcut2)
=======
        for (j = 0; (j < njcg); j++)
        {
            cg_j   = jcg[j];
            nrj    = cgindex[cg_j+1]-cgindex[cg_j];
            if (calc_image_rect(x[icg], x[cg_j], box_size, b_inv, &shift) < rcut2)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                jgid  = GET_CGINFO_GID(cginfo[cg_j]);
                if (!(i_egp_flags[jgid] & EGP_EXCL))
                {
<<<<<<< HEAD
                    add_simple(&ns_buf[jgid][shift],nrj,cg_j,
                               bHaveVdW,ngid,md,icg,jgid,cgs,bexcl,shift,fr,
=======
                    add_simple(&ns_buf[jgid][shift], nrj, cg_j,
                               bHaveVdW, ngid, md, icg, jgid, cgs, bexcl, shift, fr,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                               put_in_list);
                }
            }
        }
<<<<<<< HEAD
    } 
    else
    {
        for(j=0; (j<njcg); j++)
        {
            cg_j   = jcg[j];
            nrj    = cgindex[cg_j+1]-cgindex[cg_j];
            if ((rcut2 == 0) || (distance2(x[icg],x[cg_j]) < rcut2)) {
                jgid  = GET_CGINFO_GID(cginfo[cg_j]);
                if (!(i_egp_flags[jgid] & EGP_EXCL))
                {
                    add_simple(&ns_buf[jgid][CENTRAL],nrj,cg_j,
                               bHaveVdW,ngid,md,icg,jgid,cgs,bexcl,CENTRAL,fr,
=======
    }
    else
    {
        for (j = 0; (j < njcg); j++)
        {
            cg_j   = jcg[j];
            nrj    = cgindex[cg_j+1]-cgindex[cg_j];
            if ((rcut2 == 0) || (distance2(x[icg], x[cg_j]) < rcut2))
            {
                jgid  = GET_CGINFO_GID(cginfo[cg_j]);
                if (!(i_egp_flags[jgid] & EGP_EXCL))
                {
                    add_simple(&ns_buf[jgid][CENTRAL], nrj, cg_j,
                               bHaveVdW, ngid, md, icg, jgid, cgs, bexcl, CENTRAL, fr,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                               put_in_list);
                }
            }
        }
    }
}

/* ns_simple_core needs to be adapted for QMMM still 2005 */

static int ns_simple_core(t_forcerec *fr,
                          gmx_localtop_t *top,
                          t_mdatoms *md,
<<<<<<< HEAD
                          matrix box,rvec box_size,
                          t_excl bexcl[],atom_id *aaj,
                          int ngid,t_ns_buf **ns_buf,
                          put_in_list_t *put_in_list,gmx_bool bHaveVdW[])
{
    int      naaj,k;
    real     rlist2;
    int      nsearch,icg,jcg,igid,i0,nri,nn;
    int      *cginfo;
    t_ns_buf *nsbuf;
    /* atom_id  *i_atoms; */
    t_block  *cgs=&(top->cgs);
    t_blocka *excl=&(top->excls);
    rvec     b_inv;
    int      m;
    gmx_bool     bBox,bTriclinic;
    int      *i_egp_flags;
    
    rlist2 = sqr(fr->rlist);
    
    bBox = (fr->ePBC != epbcNONE);
    if (bBox)
    {
        for(m=0; (m<DIM); m++)
        {
            b_inv[m] = divide(1.0,box_size[m]);
=======
                          matrix box, rvec box_size,
                          t_excl bexcl[], atom_id *aaj,
                          int ngid, t_ns_buf **ns_buf,
                          put_in_list_t *put_in_list, gmx_bool bHaveVdW[])
{
    int          naaj, k;
    real         rlist2;
    int          nsearch, icg, jcg, igid, i0, nri, nn;
    int         *cginfo;
    t_ns_buf    *nsbuf;
    /* atom_id  *i_atoms; */
    t_block     *cgs  = &(top->cgs);
    t_blocka    *excl = &(top->excls);
    rvec         b_inv;
    int          m;
    gmx_bool     bBox, bTriclinic;
    int         *i_egp_flags;

    rlist2 = sqr(fr->rlist);

    bBox = (fr->ePBC != epbcNONE);
    if (bBox)
    {
        for (m = 0; (m < DIM); m++)
        {
            b_inv[m] = divide(1.0, box_size[m]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        bTriclinic = TRICLINIC(box);
    }
    else
    {
        bTriclinic = FALSE;
    }
<<<<<<< HEAD
    
    cginfo = fr->cginfo;
    
    nsearch=0;
    for (icg=fr->cg0; (icg<fr->hcg); icg++)
    {
        /*
          i0        = cgs->index[icg];
          nri       = cgs->index[icg+1]-i0;
          i_atoms   = &(cgs->a[i0]);
          i_eg_excl = fr->eg_excl + ngid*md->cENER[*i_atoms];
          setexcl(nri,i_atoms,excl,TRUE,bexcl);
        */
        igid = GET_CGINFO_GID(cginfo[icg]);
        i_egp_flags = fr->egp_flags + ngid*igid;
        setexcl(cgs->index[icg],cgs->index[icg+1],excl,TRUE,bexcl);
        
        naaj=calc_naaj(icg,cgs->nr);
        if (bTriclinic)
        {
            ns_inner_tric(fr->cg_cm,icg,i_egp_flags,naaj,&(aaj[icg]),
                          box,b_inv,rlist2,cgs,ns_buf,
                          bHaveVdW,ngid,md,bexcl,fr,put_in_list);
        }
        else
        {
            ns_inner_rect(fr->cg_cm,icg,i_egp_flags,naaj,&(aaj[icg]),
                          bBox,box_size,b_inv,rlist2,cgs,ns_buf,
                          bHaveVdW,ngid,md,bexcl,fr,put_in_list);
        }
        nsearch += naaj;
        
        for(nn=0; (nn<ngid); nn++)
        {
            for(k=0; (k<SHIFTS); k++)
=======

    cginfo = fr->cginfo;

    nsearch = 0;
    for (icg = fr->cg0; (icg < fr->hcg); icg++)
    {
        /*
           i0        = cgs->index[icg];
           nri       = cgs->index[icg+1]-i0;
           i_atoms   = &(cgs->a[i0]);
           i_eg_excl = fr->eg_excl + ngid*md->cENER[*i_atoms];
           setexcl(nri,i_atoms,excl,TRUE,bexcl);
         */
        igid        = GET_CGINFO_GID(cginfo[icg]);
        i_egp_flags = fr->egp_flags + ngid*igid;
        setexcl(cgs->index[icg], cgs->index[icg+1], excl, TRUE, bexcl);

        naaj = calc_naaj(icg, cgs->nr);
        if (bTriclinic)
        {
            ns_inner_tric(fr->cg_cm, icg, i_egp_flags, naaj, &(aaj[icg]),
                          box, b_inv, rlist2, cgs, ns_buf,
                          bHaveVdW, ngid, md, bexcl, fr, put_in_list);
        }
        else
        {
            ns_inner_rect(fr->cg_cm, icg, i_egp_flags, naaj, &(aaj[icg]),
                          bBox, box_size, b_inv, rlist2, cgs, ns_buf,
                          bHaveVdW, ngid, md, bexcl, fr, put_in_list);
        }
        nsearch += naaj;

        for (nn = 0; (nn < ngid); nn++)
        {
            for (k = 0; (k < SHIFTS); k++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                nsbuf = &(ns_buf[nn][k]);
                if (nsbuf->ncg > 0)
                {
<<<<<<< HEAD
                    put_in_list(bHaveVdW,ngid,md,icg,nn,nsbuf->ncg,nsbuf->jcg,
                                cgs->index,bexcl,k,fr,FALSE,TRUE,TRUE);
                    nsbuf->ncg=nsbuf->nj=0;
=======
                    put_in_list(bHaveVdW, ngid, md, icg, nn, nsbuf->ncg, nsbuf->jcg,
                                cgs->index, bexcl, k, fr, FALSE, TRUE, TRUE, fr->solvent_opt);
                    nsbuf->ncg = nsbuf->nj = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
            }
        }
        /* setexcl(nri,i_atoms,excl,FALSE,bexcl); */
<<<<<<< HEAD
        setexcl(cgs->index[icg],cgs->index[icg+1],excl,FALSE,bexcl);
    }
    close_neighbor_list(fr,FALSE,-1,-1,FALSE);
    
=======
        setexcl(cgs->index[icg], cgs->index[icg+1], excl, FALSE, bexcl);
    }
    close_neighbor_lists(fr, FALSE);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return nsearch;
}

/************************************************
 *
 *    N S 5     G R I D     S T U F F
 *
 ************************************************/

<<<<<<< HEAD
static inline void get_dx(int Nx,real gridx,real rc2,int xgi,real x,
                          int *dx0,int *dx1,real *dcx2)
{
    real dcx,tmp;
    int  xgi0,xgi1,i;
    
=======
static inline void get_dx(int Nx, real gridx, real rc2, int xgi, real x,
                          int *dx0, int *dx1, real *dcx2)
{
    real dcx, tmp;
    int  xgi0, xgi1, i;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (xgi < 0)
    {
        *dx0 = 0;
        xgi0 = -1;
        *dx1 = -1;
        xgi1 = 0;
    }
    else if (xgi >= Nx)
    {
        *dx0 = Nx;
        xgi0 = Nx-1;
        *dx1 = Nx-1;
        xgi1 = Nx;
    }
    else
    {
        dcx2[xgi] = 0;
<<<<<<< HEAD
        *dx0 = xgi;
        xgi0 = xgi-1;
        *dx1 = xgi;
        xgi1 = xgi+1;
    }
    
    for(i=xgi0; i>=0; i--)
=======
        *dx0      = xgi;
        xgi0      = xgi-1;
        *dx1      = xgi;
        xgi1      = xgi+1;
    }

    for (i = xgi0; i >= 0; i--)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        dcx = (i+1)*gridx-x;
        tmp = dcx*dcx;
        if (tmp >= rc2)
<<<<<<< HEAD
            break;
        *dx0 = i;
        dcx2[i] = tmp;
    }
    for(i=xgi1; i<Nx; i++)
=======
        {
            break;
        }
        *dx0    = i;
        dcx2[i] = tmp;
    }
    for (i = xgi1; i < Nx; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        dcx = i*gridx-x;
        tmp = dcx*dcx;
        if (tmp >= rc2)
        {
            break;
        }
<<<<<<< HEAD
        *dx1 = i;
=======
        *dx1    = i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        dcx2[i] = tmp;
    }
}

<<<<<<< HEAD
static inline void get_dx_dd(int Nx,real gridx,real rc2,int xgi,real x,
                             int ncpddc,int shift_min,int shift_max,
                             int *g0,int *g1,real *dcx2)
{
    real dcx,tmp;
    int  g_min,g_max,shift_home;
    
=======
static inline void get_dx_dd(int Nx, real gridx, real rc2, int xgi, real x,
                             int ncpddc, int shift_min, int shift_max,
                             int *g0, int *g1, real *dcx2)
{
    real dcx, tmp;
    int  g_min, g_max, shift_home;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (xgi < 0)
    {
        g_min = 0;
        g_max = Nx - 1;
        *g0   = 0;
        *g1   = -1;
    }
    else if (xgi >= Nx)
    {
        g_min = 0;
        g_max = Nx - 1;
        *g0   = Nx;
        *g1   = Nx - 1;
    }
    else
    {
        if (ncpddc == 0)
        {
            g_min = 0;
            g_max = Nx - 1;
        }
        else
        {
            if (xgi < ncpddc)
            {
                shift_home = 0;
            }
            else
            {
                shift_home = -1;
            }
            g_min = (shift_min == shift_home ? 0          : ncpddc);
            g_max = (shift_max == shift_home ? ncpddc - 1 : Nx - 1);
        }
        if (shift_min > 0)
        {
            *g0 = g_min;
            *g1 = g_min - 1;
        }
        else if (shift_max < 0)
        {
            *g0 = g_max + 1;
            *g1 = g_max;
        }
        else
        {
<<<<<<< HEAD
            *g0 = xgi;
            *g1 = xgi;
            dcx2[xgi] = 0;
        }
    }
    
=======
            *g0       = xgi;
            *g1       = xgi;
            dcx2[xgi] = 0;
        }
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    while (*g0 > g_min)
    {
        /* Check one grid cell down */
        dcx = ((*g0 - 1) + 1)*gridx - x;
        tmp = dcx*dcx;
        if (tmp >= rc2)
        {
            break;
        }
        (*g0)--;
        dcx2[*g0] = tmp;
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    while (*g1 < g_max)
    {
        /* Check one grid cell up */
        dcx = (*g1 + 1)*gridx - x;
        tmp = dcx*dcx;
        if (tmp >= rc2)
        {
            break;
        }
        (*g1)++;
        dcx2[*g1] = tmp;
    }
}


#define sqr(x) ((x)*(x))
<<<<<<< HEAD
#define calc_dx2(XI,YI,ZI,y) (sqr(XI-y[XX]) + sqr(YI-y[YY]) + sqr(ZI-y[ZZ]))
#define calc_cyl_dx2(XI,YI,y) (sqr(XI-y[XX]) + sqr(YI-y[YY]))
=======
#define calc_dx2(XI, YI, ZI, y) (sqr(XI-y[XX]) + sqr(YI-y[YY]) + sqr(ZI-y[ZZ]))
#define calc_cyl_dx2(XI, YI, y) (sqr(XI-y[XX]) + sqr(YI-y[YY]))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/****************************************************
 *
 *    F A S T   N E I G H B O R  S E A R C H I N G
 *
<<<<<<< HEAD
 *    Optimized neighboursearching routine using grid 
=======
 *    Optimized neighboursearching routine using grid
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *    at least 1x1x1, see GROMACS manual
 *
 ****************************************************/

<<<<<<< HEAD
static void do_longrange(t_commrec *cr,gmx_localtop_t *top,t_forcerec *fr,
                         int ngid,t_mdatoms *md,int icg,
                         int jgid,int nlr,
                         atom_id lr[],t_excl bexcl[],int shift,
                         rvec x[],rvec box_size,t_nrnb *nrnb,
                         real lambda,real *dvdlambda,
                         gmx_grppairener_t *grppener,
                         gmx_bool bDoVdW,gmx_bool bDoCoul,
                         gmx_bool bEvaluateNow,put_in_list_t *put_in_list,
                         gmx_bool bHaveVdW[],
                         gmx_bool bDoForces,rvec *f)
{
    int n,i;
    t_nblist *nl;
    
    for(n=0; n<fr->nnblists; n++)
    {
        for(i=0; (i<eNL_NR); i++)
        {
            nl = &fr->nblists[n].nlist_lr[i];
            if ((nl->nri > nl->maxnri-32) || bEvaluateNow)
            {
                close_neighbor_list(fr,TRUE,n,i,FALSE);
                /* Evaluate the energies and forces */
                do_nonbonded(cr,fr,x,f,md,NULL,
                             grppener->ener[fr->bBHAM ? egBHAMLR : egLJLR],
                             grppener->ener[egCOULLR],
							 grppener->ener[egGB],box_size,
                             nrnb,lambda,dvdlambda,n,i,
                             GMX_DONB_LR | GMX_DONB_FORCES);
                
                reset_neighbor_list(fr,TRUE,n,i);
            }
        }
    }
    
    if (!bEvaluateNow)
    {  
        /* Put the long range particles in a list */
        /* do_longrange is never called for QMMM  */
        put_in_list(bHaveVdW,ngid,md,icg,jgid,nlr,lr,top->cgs.index,
                    bexcl,shift,fr,TRUE,bDoVdW,bDoCoul);
    }
}

static void get_cutoff2(t_forcerec *fr,gmx_bool bDoLongRange,
                        real *rvdw2,real *rcoul2,
                        real *rs2,real *rm2,real *rl2)
{
    *rs2 = sqr(fr->rlist);
=======

static void get_cutoff2(t_forcerec *fr, gmx_bool bDoLongRange,
                        real *rvdw2, real *rcoul2,
                        real *rs2, real *rm2, real *rl2)
{
    *rs2 = sqr(fr->rlist);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (bDoLongRange && fr->bTwinRange)
    {
        /* The VdW and elec. LR cut-off's could be different,
         * so we can not simply set them to rlistlong.
         */
        if (EVDW_MIGHT_BE_ZERO_AT_CUTOFF(fr->vdwtype) &&
            fr->rvdw > fr->rlist)
        {
            *rvdw2  = sqr(fr->rlistlong);
        }
        else
        {
            *rvdw2  = sqr(fr->rvdw);
        }
        if (EEL_MIGHT_BE_ZERO_AT_CUTOFF(fr->eeltype) &&
            fr->rcoulomb > fr->rlist)
        {
            *rcoul2 = sqr(fr->rlistlong);
        }
        else
        {
            *rcoul2 = sqr(fr->rcoulomb);
        }
    }
    else
    {
        /* Workaround for a gcc -O3 or -ffast-math problem */
        *rvdw2  = *rs2;
        *rcoul2 = *rs2;
    }
<<<<<<< HEAD
    *rm2 = min(*rvdw2,*rcoul2);
    *rl2 = max(*rvdw2,*rcoul2);
}

static void init_nsgrid_lists(t_forcerec *fr,int ngid,gmx_ns_t *ns)
{
    real rvdw2,rcoul2,rs2,rm2,rl2;
    int j;

    get_cutoff2(fr,TRUE,&rvdw2,&rcoul2,&rs2,&rm2,&rl2);

    /* Short range buffers */
    snew(ns->nl_sr,ngid);
    /* Counters */
    snew(ns->nsr,ngid);
    snew(ns->nlr_ljc,ngid);
    snew(ns->nlr_one,ngid);
    
    if (rm2 > rs2)
    {
            /* Long range VdW and Coul buffers */
        snew(ns->nl_lr_ljc,ngid);
    }
    if (rl2 > rm2)
    {
        /* Long range VdW or Coul only buffers */
        snew(ns->nl_lr_one,ngid);
    }
    for(j=0; (j<ngid); j++) {
        snew(ns->nl_sr[j],MAX_CG);
        if (rm2 > rs2)
        {
            snew(ns->nl_lr_ljc[j],MAX_CG);
        }
        if (rl2 > rm2)
        {
            snew(ns->nl_lr_one[j],MAX_CG);
        }
=======
    *rm2 = min(*rvdw2, *rcoul2);
    *rl2 = max(*rvdw2, *rcoul2);
}

static void init_nsgrid_lists(t_forcerec *fr, int ngid, gmx_ns_t *ns)
{
    real rvdw2, rcoul2, rs2, rm2, rl2;
    int  j;

    get_cutoff2(fr, TRUE, &rvdw2, &rcoul2, &rs2, &rm2, &rl2);

    /* Short range buffers */
    snew(ns->nl_sr, ngid);
    /* Counters */
    snew(ns->nsr, ngid);
    snew(ns->nlr_ljc, ngid);
    snew(ns->nlr_one, ngid);

    /* Always allocate both list types, since rcoulomb might now change with PME load balancing */
    /* Long range VdW and Coul buffers */
    snew(ns->nl_lr_ljc, ngid);
    /* Long range VdW or Coul only buffers */
    snew(ns->nl_lr_one, ngid);

    for (j = 0; (j < ngid); j++)
    {
        snew(ns->nl_sr[j], MAX_CG);
        snew(ns->nl_lr_ljc[j], MAX_CG);
        snew(ns->nl_lr_one[j], MAX_CG);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    if (debug)
    {
        fprintf(debug,
                "ns5_core: rs2 = %g, rm2 = %g, rl2 = %g (nm^2)\n",
<<<<<<< HEAD
                rs2,rm2,rl2);
    }
}

static int nsgrid_core(FILE *log,t_commrec *cr,t_forcerec *fr,
                       matrix box,rvec box_size,int ngid,
                       gmx_localtop_t *top,
                       t_grid *grid,rvec x[],
                       t_excl bexcl[],gmx_bool *bExcludeAlleg,
                       t_nrnb *nrnb,t_mdatoms *md,
                       real lambda,real *dvdlambda,
                       gmx_grppairener_t *grppener,
                       put_in_list_t *put_in_list,
                       gmx_bool bHaveVdW[],
                       gmx_bool bDoLongRange,gmx_bool bDoForces,rvec *f,
                       gmx_bool bMakeQMMMnblist)
{
    gmx_ns_t *ns;
    atom_id **nl_lr_ljc,**nl_lr_one,**nl_sr;
    int     *nlr_ljc,*nlr_one,*nsr;
    gmx_domdec_t *dd=NULL;
    t_block *cgs=&(top->cgs);
    int     *cginfo=fr->cginfo;
    /* atom_id *i_atoms,*cgsindex=cgs->index; */
    ivec    sh0,sh1,shp;
    int     cell_x,cell_y,cell_z;
    int     d,tx,ty,tz,dx,dy,dz,cj;
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
    int     zsh_ty,zsh_tx,ysh_tx;
#endif
    int     dx0,dx1,dy0,dy1,dz0,dz1;
    int     Nx,Ny,Nz,shift=-1,j,nrj,nns,nn=-1;
    real    gridx,gridy,gridz,grid_x,grid_y,grid_z;
    real    *dcx2,*dcy2,*dcz2;
    int     zgi,ygi,xgi;
    int     cg0,cg1,icg=-1,cgsnr,i0,igid,nri,naaj,max_jcg;
    int     jcg0,jcg1,jjcg,cgj0,jgid;
    int     *grida,*gridnra,*gridind;
    gmx_bool    rvdw_lt_rcoul,rcoul_lt_rvdw;
    rvec    xi,*cgcm,grid_offset;
    real    r2,rs2,rvdw2,rcoul2,rm2,rl2,XI,YI,ZI,dcx,dcy,dcz,tmp1,tmp2;
    int     *i_egp_flags;
    gmx_bool    bDomDec,bTriclinicX,bTriclinicY;
    ivec    ncpddc;
    
    ns = &fr->ns;
    
=======
                rs2, rm2, rl2);
    }
}

static int nsgrid_core(FILE *log, t_commrec *cr, t_forcerec *fr,
                       matrix box, rvec box_size, int ngid,
                       gmx_localtop_t *top,
                       t_grid *grid, rvec x[],
                       t_excl bexcl[], gmx_bool *bExcludeAlleg,
                       t_nrnb *nrnb, t_mdatoms *md,
                       real *lambda, real *dvdlambda,
                       gmx_grppairener_t *grppener,
                       put_in_list_t *put_in_list,
                       gmx_bool bHaveVdW[],
                       gmx_bool bDoLongRange, gmx_bool bMakeQMMMnblist)
{
    gmx_ns_t     *ns;
    atom_id     **nl_lr_ljc, **nl_lr_one, **nl_sr;
    int          *nlr_ljc, *nlr_one, *nsr;
    gmx_domdec_t *dd     = NULL;
    t_block      *cgs    = &(top->cgs);
    int          *cginfo = fr->cginfo;
    /* atom_id *i_atoms,*cgsindex=cgs->index; */
    ivec          sh0, sh1, shp;
    int           cell_x, cell_y, cell_z;
    int           d, tx, ty, tz, dx, dy, dz, cj;
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
    int           zsh_ty, zsh_tx, ysh_tx;
#endif
    int           dx0, dx1, dy0, dy1, dz0, dz1;
    int           Nx, Ny, Nz, shift = -1, j, nrj, nns, nn = -1;
    real          gridx, gridy, gridz, grid_x, grid_y, grid_z;
    real         *dcx2, *dcy2, *dcz2;
    int           zgi, ygi, xgi;
    int           cg0, cg1, icg = -1, cgsnr, i0, igid, nri, naaj, max_jcg;
    int           jcg0, jcg1, jjcg, cgj0, jgid;
    int          *grida, *gridnra, *gridind;
    gmx_bool      rvdw_lt_rcoul, rcoul_lt_rvdw;
    rvec          xi, *cgcm, grid_offset;
    real          r2, rs2, rvdw2, rcoul2, rm2, rl2, XI, YI, ZI, dcx, dcy, dcz, tmp1, tmp2;
    int          *i_egp_flags;
    gmx_bool      bDomDec, bTriclinicX, bTriclinicY;
    ivec          ncpddc;

    ns = &fr->ns;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    bDomDec = DOMAINDECOMP(cr);
    if (bDomDec)
    {
        dd = cr->dd;
    }
<<<<<<< HEAD
    
    bTriclinicX = ((YY < grid->npbcdim &&
                    (!bDomDec || dd->nc[YY]==1) && box[YY][XX] != 0) ||
                   (ZZ < grid->npbcdim &&
                    (!bDomDec || dd->nc[ZZ]==1) && box[ZZ][XX] != 0));
    bTriclinicY =  (ZZ < grid->npbcdim &&
                    (!bDomDec || dd->nc[ZZ]==1) && box[ZZ][YY] != 0);
    
    cgsnr    = cgs->nr;

    get_cutoff2(fr,bDoLongRange,&rvdw2,&rcoul2,&rs2,&rm2,&rl2);

    rvdw_lt_rcoul = (rvdw2 >= rcoul2);
    rcoul_lt_rvdw = (rcoul2 >= rvdw2);
    
=======

    bTriclinicX = ((YY < grid->npbcdim &&
                    (!bDomDec || dd->nc[YY] == 1) && box[YY][XX] != 0) ||
                   (ZZ < grid->npbcdim &&
                    (!bDomDec || dd->nc[ZZ] == 1) && box[ZZ][XX] != 0));
    bTriclinicY =  (ZZ < grid->npbcdim &&
                    (!bDomDec || dd->nc[ZZ] == 1) && box[ZZ][YY] != 0);

    cgsnr    = cgs->nr;

    get_cutoff2(fr, bDoLongRange, &rvdw2, &rcoul2, &rs2, &rm2, &rl2);

    rvdw_lt_rcoul = (rvdw2 >= rcoul2);
    rcoul_lt_rvdw = (rcoul2 >= rvdw2);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (bMakeQMMMnblist)
    {
        rm2 = rl2;
        rs2 = rl2;
    }

    nl_sr     = ns->nl_sr;
    nsr       = ns->nsr;
    nl_lr_ljc = ns->nl_lr_ljc;
    nl_lr_one = ns->nl_lr_one;
    nlr_ljc   = ns->nlr_ljc;
    nlr_one   = ns->nlr_one;
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Unpack arrays */
    cgcm    = fr->cg_cm;
    Nx      = grid->n[XX];
    Ny      = grid->n[YY];
    Nz      = grid->n[ZZ];
    grida   = grid->a;
    gridind = grid->index;
    gridnra = grid->nra;
    nns     = 0;
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gridx      = grid->cell_size[XX];
    gridy      = grid->cell_size[YY];
    gridz      = grid->cell_size[ZZ];
    grid_x     = 1/gridx;
    grid_y     = 1/gridy;
    grid_z     = 1/gridz;
<<<<<<< HEAD
    copy_rvec(grid->cell_offset,grid_offset);
    copy_ivec(grid->ncpddc,ncpddc);
    dcx2       = grid->dcx2;
    dcy2       = grid->dcy2;
    dcz2       = grid->dcz2;
    
=======
    copy_rvec(grid->cell_offset, grid_offset);
    copy_ivec(grid->ncpddc, ncpddc);
    dcx2       = grid->dcx2;
    dcy2       = grid->dcy2;
    dcz2       = grid->dcz2;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
    zsh_ty = floor(-box[ZZ][YY]/box[YY][YY]+0.5);
    zsh_tx = floor(-box[ZZ][XX]/box[XX][XX]+0.5);
    ysh_tx = floor(-box[YY][XX]/box[XX][XX]+0.5);
<<<<<<< HEAD
    if (zsh_tx!=0 && ysh_tx!=0)
=======
    if (zsh_tx != 0 && ysh_tx != 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* This could happen due to rounding, when both ratios are 0.5 */
        ysh_tx = 0;
    }
#endif
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    debug_gmx();

    if (fr->n_tpi)
    {
        /* We only want a list for the test particle */
        cg0 = cgsnr - 1;
    }
    else
    {
        cg0 = grid->icg0;
    }
    cg1 = grid->icg1;

    /* Set the shift range */
<<<<<<< HEAD
    for(d=0; d<DIM; d++)
=======
    for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        sh0[d] = -1;
        sh1[d] = 1;
        /* Check if we need periodicity shifts.
         * Without PBC or with domain decomposition we don't need them.
         */
        if (d >= ePBC2npbcdim(fr->ePBC) || (bDomDec && dd->nc[d] > 1))
        {
            shp[d] = 0;
        }
        else
        {
            if (d == XX &&
                box[XX][XX] - fabs(box[YY][XX]) - fabs(box[ZZ][XX]) < sqrt(rl2))
            {
                shp[d] = 2;
            }
            else
            {
                shp[d] = 1;
            }
        }
    }
<<<<<<< HEAD
    
    /* Loop over charge groups */
    for(icg=cg0; (icg < cg1); icg++)
=======

    /* Loop over charge groups */
    for (icg = cg0; (icg < cg1); icg++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        igid = GET_CGINFO_GID(cginfo[icg]);
        /* Skip this charge group if all energy groups are excluded! */
        if (bExcludeAlleg[igid])
        {
            continue;
        }
<<<<<<< HEAD
        
        i0   = cgs->index[icg];
        
        if (bMakeQMMMnblist)
        { 
            /* Skip this charge group if it is not a QM atom while making a
             * QM/MM neighbourlist
             */
            if (md->bQM[i0]==FALSE)
            {
                continue; /* MM particle, go to next particle */ 
            }
            
            /* Compute the number of charge groups that fall within the control
             * of this one (icg)
             */
            naaj    = calc_naaj(icg,cgsnr);
            jcg0    = icg;
            jcg1    = icg + naaj;
            max_jcg = cgsnr;       
        } 
        else
        { 
            /* make a normal neighbourlist */
            
            if (bDomDec)
            {
                /* Get the j charge-group and dd cell shift ranges */
                dd_get_ns_ranges(cr->dd,icg,&jcg0,&jcg1,sh0,sh1);
=======

        i0   = cgs->index[icg];

        if (bMakeQMMMnblist)
        {
            /* Skip this charge group if it is not a QM atom while making a
             * QM/MM neighbourlist
             */
            if (md->bQM[i0] == FALSE)
            {
                continue; /* MM particle, go to next particle */
            }

            /* Compute the number of charge groups that fall within the control
             * of this one (icg)
             */
            naaj    = calc_naaj(icg, cgsnr);
            jcg0    = icg;
            jcg1    = icg + naaj;
            max_jcg = cgsnr;
        }
        else
        {
            /* make a normal neighbourlist */

            if (bDomDec)
            {
                /* Get the j charge-group and dd cell shift ranges */
                dd_get_ns_ranges(cr->dd, icg, &jcg0, &jcg1, sh0, sh1);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                max_jcg = 0;
            }
            else
            {
                /* Compute the number of charge groups that fall within the control
                 * of this one (icg)
                 */
<<<<<<< HEAD
                naaj = calc_naaj(icg,cgsnr);
                jcg0 = icg;
                jcg1 = icg + naaj;
                
=======
                naaj = calc_naaj(icg, cgsnr);
                jcg0 = icg;
                jcg1 = icg + naaj;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                if (fr->n_tpi)
                {
                    /* The i-particle is awlways the test particle,
                     * so we want all j-particles
                     */
                    max_jcg = cgsnr - 1;
                }
                else
                {
                    max_jcg  = jcg1 - cgsnr;
                }
            }
        }
<<<<<<< HEAD
        
        i_egp_flags = fr->egp_flags + igid*ngid;
        
        /* Set the exclusions for the atoms in charge group icg using a bitmask */
        setexcl(i0,cgs->index[icg+1],&top->excls,TRUE,bexcl);
        
        ci2xyz(grid,icg,&cell_x,&cell_y,&cell_z);
        
        /* Changed iicg to icg, DvdS 990115 
         * (but see consistency check above, DvdS 990330) 
         */
#ifdef NS5DB
        fprintf(log,"icg=%5d, naaj=%5d, cell %d %d %d\n",
                icg,naaj,cell_x,cell_y,cell_z);
#endif
        /* Loop over shift vectors in three dimensions */
        for (tz=-shp[ZZ]; tz<=shp[ZZ]; tz++)
=======

        i_egp_flags = fr->egp_flags + igid*ngid;

        /* Set the exclusions for the atoms in charge group icg using a bitmask */
        setexcl(i0, cgs->index[icg+1], &top->excls, TRUE, bexcl);

        ci2xyz(grid, icg, &cell_x, &cell_y, &cell_z);

        /* Changed iicg to icg, DvdS 990115
         * (but see consistency check above, DvdS 990330)
         */
#ifdef NS5DB
        fprintf(log, "icg=%5d, naaj=%5d, cell %d %d %d\n",
                icg, naaj, cell_x, cell_y, cell_z);
#endif
        /* Loop over shift vectors in three dimensions */
        for (tz = -shp[ZZ]; tz <= shp[ZZ]; tz++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            ZI = cgcm[icg][ZZ]+tz*box[ZZ][ZZ];
            /* Calculate range of cells in Z direction that have the shift tz */
            zgi = cell_z + tz*Nz;
#define FAST_DD_NS
#ifndef FAST_DD_NS
<<<<<<< HEAD
            get_dx(Nz,gridz,rl2,zgi,ZI,&dz0,&dz1,dcz2);
#else
            get_dx_dd(Nz,gridz,rl2,zgi,ZI-grid_offset[ZZ],
                      ncpddc[ZZ],sh0[ZZ],sh1[ZZ],&dz0,&dz1,dcz2);
=======
            get_dx(Nz, gridz, rl2, zgi, ZI, &dz0, &dz1, dcz2);
#else
            get_dx_dd(Nz, gridz, rl2, zgi, ZI-grid_offset[ZZ],
                      ncpddc[ZZ], sh0[ZZ], sh1[ZZ], &dz0, &dz1, dcz2);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
            if (dz0 > dz1)
            {
                continue;
            }
<<<<<<< HEAD
            for (ty=-shp[YY]; ty<=shp[YY]; ty++)
=======
            for (ty = -shp[YY]; ty <= shp[YY]; ty++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                YI = cgcm[icg][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
                /* Calculate range of cells in Y direction that have the shift ty */
                if (bTriclinicY)
                {
                    ygi = (int)(Ny + (YI - grid_offset[YY])*grid_y) - Ny;
                }
                else
                {
                    ygi = cell_y + ty*Ny;
                }
#ifndef FAST_DD_NS
<<<<<<< HEAD
                get_dx(Ny,gridy,rl2,ygi,YI,&dy0,&dy1,dcy2);
#else
                get_dx_dd(Ny,gridy,rl2,ygi,YI-grid_offset[YY],
                          ncpddc[YY],sh0[YY],sh1[YY],&dy0,&dy1,dcy2);
=======
                get_dx(Ny, gridy, rl2, ygi, YI, &dy0, &dy1, dcy2);
#else
                get_dx_dd(Ny, gridy, rl2, ygi, YI-grid_offset[YY],
                          ncpddc[YY], sh0[YY], sh1[YY], &dy0, &dy1, dcy2);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
                if (dy0 > dy1)
                {
                    continue;
                }
<<<<<<< HEAD
                for (tx=-shp[XX]; tx<=shp[XX]; tx++)
=======
                for (tx = -shp[XX]; tx <= shp[XX]; tx++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    XI = cgcm[icg][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
                    /* Calculate range of cells in X direction that have the shift tx */
                    if (bTriclinicX)
                    {
                        xgi = (int)(Nx + (XI - grid_offset[XX])*grid_x) - Nx;
                    }
                    else
                    {
                        xgi = cell_x + tx*Nx;
                    }
#ifndef FAST_DD_NS
<<<<<<< HEAD
                    get_dx(Nx,gridx,rl2,xgi*Nx,XI,&dx0,&dx1,dcx2);
#else
                    get_dx_dd(Nx,gridx,rl2,xgi,XI-grid_offset[XX],
                              ncpddc[XX],sh0[XX],sh1[XX],&dx0,&dx1,dcx2);
=======
                    get_dx(Nx, gridx, rl2, xgi*Nx, XI, &dx0, &dx1, dcx2);
#else
                    get_dx_dd(Nx, gridx, rl2, xgi, XI-grid_offset[XX],
                              ncpddc[XX], sh0[XX], sh1[XX], &dx0, &dx1, dcx2);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
                    if (dx0 > dx1)
                    {
                        continue;
                    }
                    /* Adress: an explicit cg that has a weigthing function of 0 is excluded
                     *  from the neigbour list as it will not interact  */
<<<<<<< HEAD
                    /*if (fr->adress_type != eAdressOff){
                        if (md->wf[cgs->index[icg]]==0 && egp_explicit(fr, igid)){
                            continue;
                        }
                    }*/
                    /* Get shift vector */	  
                    shift=XYZ2IS(tx,ty,tz);
#ifdef NS5DB
                    range_check(shift,0,SHIFTS);
#endif
                    for(nn=0; (nn<ngid); nn++)
                    {
                        nsr[nn]      = 0;
                        nlr_ljc[nn]  = 0;
                        nlr_one[nn] = 0;
                    }
#ifdef NS5DB
                    fprintf(log,"shift: %2d, dx0,1: %2d,%2d, dy0,1: %2d,%2d, dz0,1: %2d,%2d\n",
                            shift,dx0,dx1,dy0,dy1,dz0,dz1);
                    fprintf(log,"cgcm: %8.3f  %8.3f  %8.3f\n",cgcm[icg][XX],
                            cgcm[icg][YY],cgcm[icg][ZZ]);
                    fprintf(log,"xi:   %8.3f  %8.3f  %8.3f\n",XI,YI,ZI);
#endif
                    for (dx=dx0; (dx<=dx1); dx++)
                    {
                        tmp1 = rl2 - dcx2[dx];
                        for (dy=dy0; (dy<=dy1); dy++)
=======
                    if (fr->adress_type != eAdressOff)
                    {
                        if (md->wf[cgs->index[icg]] <= GMX_REAL_EPS && egp_explicit(fr, igid))
                        {
                            continue;
                        }
                    }
                    /* Get shift vector */
                    shift = XYZ2IS(tx, ty, tz);
#ifdef NS5DB
                    range_check(shift, 0, SHIFTS);
#endif
                    for (nn = 0; (nn < ngid); nn++)
                    {
                        nsr[nn]      = 0;
                        nlr_ljc[nn]  = 0;
                        nlr_one[nn]  = 0;
                    }
#ifdef NS5DB
                    fprintf(log, "shift: %2d, dx0,1: %2d,%2d, dy0,1: %2d,%2d, dz0,1: %2d,%2d\n",
                            shift, dx0, dx1, dy0, dy1, dz0, dz1);
                    fprintf(log, "cgcm: %8.3f  %8.3f  %8.3f\n", cgcm[icg][XX],
                            cgcm[icg][YY], cgcm[icg][ZZ]);
                    fprintf(log, "xi:   %8.3f  %8.3f  %8.3f\n", XI, YI, ZI);
#endif
                    for (dx = dx0; (dx <= dx1); dx++)
                    {
                        tmp1 = rl2 - dcx2[dx];
                        for (dy = dy0; (dy <= dy1); dy++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                        {
                            tmp2 = tmp1 - dcy2[dy];
                            if (tmp2 > 0)
                            {
<<<<<<< HEAD
                                for (dz=dz0; (dz<=dz1); dz++) {
                                    if (tmp2 > dcz2[dz]) {
                                        /* Find grid-cell cj in which possible neighbours are */
                                        cj   = xyz2ci(Ny,Nz,dx,dy,dz);
                                        
                                        /* Check out how many cgs (nrj) there in this cell */
                                        nrj  = gridnra[cj];
                                        
                                        /* Find the offset in the cg list */
                                        cgj0 = gridind[cj];
                                        
=======
                                for (dz = dz0; (dz <= dz1); dz++)
                                {
                                    if (tmp2 > dcz2[dz])
                                    {
                                        /* Find grid-cell cj in which possible neighbours are */
                                        cj   = xyz2ci(Ny, Nz, dx, dy, dz);

                                        /* Check out how many cgs (nrj) there in this cell */
                                        nrj  = gridnra[cj];

                                        /* Find the offset in the cg list */
                                        cgj0 = gridind[cj];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                        /* Check if all j's are out of range so we
                                         * can skip the whole cell.
                                         * Should save some time, especially with DD.
                                         */
                                        if (nrj == 0 ||
                                            (grida[cgj0] >= max_jcg &&
                                             (grida[cgj0] >= jcg1 || grida[cgj0+nrj-1] < jcg0)))
                                        {
                                            continue;
                                        }
<<<<<<< HEAD
                                        
                                        /* Loop over cgs */
                                        for (j=0; (j<nrj); j++)
                                        {
                                            jjcg = grida[cgj0+j];
                                            
=======

                                        /* Loop over cgs */
                                        for (j = 0; (j < nrj); j++)
                                        {
                                            jjcg = grida[cgj0+j];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                            /* check whether this guy is in range! */
                                            if ((jjcg >= jcg0 && jjcg < jcg1) ||
                                                (jjcg < max_jcg))
                                            {
<<<<<<< HEAD
                                                r2=calc_dx2(XI,YI,ZI,cgcm[jjcg]);
                                                if (r2 < rl2) {
=======
                                                r2 = calc_dx2(XI, YI, ZI, cgcm[jjcg]);
                                                if (r2 < rl2)
                                                {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                                    /* jgid = gid[cgsatoms[cgsindex[jjcg]]]; */
                                                    jgid = GET_CGINFO_GID(cginfo[jjcg]);
                                                    /* check energy group exclusions */
                                                    if (!(i_egp_flags[jgid] & EGP_EXCL))
                                                    {
                                                        if (r2 < rs2)
                                                        {
                                                            if (nsr[jgid] >= MAX_CG)
                                                            {
<<<<<<< HEAD
                                                                put_in_list(bHaveVdW,ngid,md,icg,jgid,
                                                                            nsr[jgid],nl_sr[jgid],
                                                                            cgs->index,/* cgsatoms, */ bexcl,
                                                                            shift,fr,FALSE,TRUE,TRUE);
                                                                nsr[jgid]=0;
                                                            }
                                                            nl_sr[jgid][nsr[jgid]++]=jjcg;
                                                        } 
=======
                                                                /* Add to short-range list */
                                                                put_in_list(bHaveVdW, ngid, md, icg, jgid,
                                                                            nsr[jgid], nl_sr[jgid],
                                                                            cgs->index, /* cgsatoms, */ bexcl,
                                                                            shift, fr, FALSE, TRUE, TRUE, fr->solvent_opt);
                                                                nsr[jgid] = 0;
                                                            }
                                                            nl_sr[jgid][nsr[jgid]++] = jjcg;
                                                        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                                        else if (r2 < rm2)
                                                        {
                                                            if (nlr_ljc[jgid] >= MAX_CG)
                                                            {
<<<<<<< HEAD
                                                                do_longrange(cr,top,fr,ngid,md,icg,jgid,
                                                                             nlr_ljc[jgid],
                                                                             nl_lr_ljc[jgid],bexcl,shift,x,
                                                                             box_size,nrnb,
                                                                             lambda,dvdlambda,
                                                                             grppener,
                                                                             TRUE,TRUE,FALSE,
                                                                             put_in_list,
                                                                             bHaveVdW,
                                                                             bDoForces,f);
                                                                nlr_ljc[jgid]=0;
                                                            }
                                                            nl_lr_ljc[jgid][nlr_ljc[jgid]++]=jjcg;
                                                        }
                                                        else
                                                        {
                                                            if (nlr_one[jgid] >= MAX_CG) {
                                                                do_longrange(cr,top,fr,ngid,md,icg,jgid,
                                                                             nlr_one[jgid],
                                                                             nl_lr_one[jgid],bexcl,shift,x,
                                                                             box_size,nrnb,
                                                                             lambda,dvdlambda,
                                                                             grppener,
                                                                             rvdw_lt_rcoul,rcoul_lt_rvdw,FALSE,
                                                                             put_in_list,
                                                                             bHaveVdW,
                                                                             bDoForces,f);
                                                                nlr_one[jgid]=0;
                                                            }
                                                            nl_lr_one[jgid][nlr_one[jgid]++]=jjcg;
=======
                                                                /* Add to LJ+coulomb long-range list */
                                                                put_in_list(bHaveVdW, ngid, md, icg, jgid,
                                                                            nlr_ljc[jgid], nl_lr_ljc[jgid], top->cgs.index,
                                                                            bexcl, shift, fr, TRUE, TRUE, TRUE, fr->solvent_opt);
                                                                nlr_ljc[jgid] = 0;
                                                            }
                                                            nl_lr_ljc[jgid][nlr_ljc[jgid]++] = jjcg;
                                                        }
                                                        else
                                                        {
                                                            if (nlr_one[jgid] >= MAX_CG)
                                                            {
                                                                /* Add to long-range list with only coul, or only LJ */
                                                                put_in_list(bHaveVdW, ngid, md, icg, jgid,
                                                                            nlr_one[jgid], nl_lr_one[jgid], top->cgs.index,
                                                                            bexcl, shift, fr, TRUE, rvdw_lt_rcoul, rcoul_lt_rvdw, fr->solvent_opt);
                                                                nlr_one[jgid] = 0;
                                                            }
                                                            nl_lr_one[jgid][nlr_one[jgid]++] = jjcg;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                                        }
                                                    }
                                                }
                                                nns++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* CHECK whether there is anything left in the buffers */
<<<<<<< HEAD
                    for(nn=0; (nn<ngid); nn++)
                    {
                        if (nsr[nn] > 0)
                        {
                            put_in_list(bHaveVdW,ngid,md,icg,nn,nsr[nn],nl_sr[nn],
                                        cgs->index, /* cgsatoms, */ bexcl,
                                        shift,fr,FALSE,TRUE,TRUE);
                        }
                        
                        if (nlr_ljc[nn] > 0)
                        {
                            do_longrange(cr,top,fr,ngid,md,icg,nn,nlr_ljc[nn],
                                         nl_lr_ljc[nn],bexcl,shift,x,box_size,nrnb,
                                         lambda,dvdlambda,grppener,TRUE,TRUE,FALSE,
                                         put_in_list,bHaveVdW,bDoForces,f);
                        }
                        
                        if (nlr_one[nn] > 0)
                        {
                            do_longrange(cr,top,fr,ngid,md,icg,nn,nlr_one[nn],
                                         nl_lr_one[nn],bexcl,shift,x,box_size,nrnb,
                                         lambda,dvdlambda,grppener,
                                         rvdw_lt_rcoul,rcoul_lt_rvdw,FALSE,
                                         put_in_list,bHaveVdW,bDoForces,f);
=======
                    for (nn = 0; (nn < ngid); nn++)
                    {
                        if (nsr[nn] > 0)
                        {
                            put_in_list(bHaveVdW, ngid, md, icg, nn, nsr[nn], nl_sr[nn],
                                        cgs->index, /* cgsatoms, */ bexcl,
                                        shift, fr, FALSE, TRUE, TRUE, fr->solvent_opt);
                        }

                        if (nlr_ljc[nn] > 0)
                        {
                            put_in_list(bHaveVdW, ngid, md, icg, nn, nlr_ljc[nn],
                                        nl_lr_ljc[nn], top->cgs.index,
                                        bexcl, shift, fr, TRUE, TRUE, TRUE, fr->solvent_opt);
                        }

                        if (nlr_one[nn] > 0)
                        {
                            put_in_list(bHaveVdW, ngid, md, icg, nn, nlr_one[nn],
                                        nl_lr_one[nn], top->cgs.index,
                                        bexcl, shift, fr, TRUE, rvdw_lt_rcoul, rcoul_lt_rvdw, fr->solvent_opt);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                        }
                    }
                }
            }
        }
        /* setexcl(nri,i_atoms,&top->atoms.excl,FALSE,bexcl); */
<<<<<<< HEAD
        setexcl(cgs->index[icg],cgs->index[icg+1],&top->excls,FALSE,bexcl);
    }
    /* Perform any left over force calculations */
    for (nn=0; (nn<ngid); nn++)
    {
        if (rm2 > rs2)
        {
            do_longrange(cr,top,fr,0,md,icg,nn,nlr_ljc[nn],
                         nl_lr_ljc[nn],bexcl,shift,x,box_size,nrnb,
                         lambda,dvdlambda,grppener,
                         TRUE,TRUE,TRUE,put_in_list,bHaveVdW,bDoForces,f);
        }
        if (rl2 > rm2) {
            do_longrange(cr,top,fr,0,md,icg,nn,nlr_one[nn],
                         nl_lr_one[nn],bexcl,shift,x,box_size,nrnb,
                         lambda,dvdlambda,grppener,
                         rvdw_lt_rcoul,rcoul_lt_rvdw,
                         TRUE,put_in_list,bHaveVdW,bDoForces,f);
        }
    }
    debug_gmx();
    
    /* Close off short range neighbourlists */
    close_neighbor_list(fr,FALSE,-1,-1,bMakeQMMMnblist);
    
    return nns;
}

void ns_realloc_natoms(gmx_ns_t *ns,int natoms)
{
    int i;
    
    if (natoms > ns->nra_alloc)
    {
        ns->nra_alloc = over_alloc_dd(natoms);
        srenew(ns->bexcl,ns->nra_alloc);
        for(i=0; i<ns->nra_alloc; i++)
=======
        setexcl(cgs->index[icg], cgs->index[icg+1], &top->excls, FALSE, bexcl);
    }
    /* No need to perform any left-over force calculations anymore (as we used to do here)
     * since we now save the proper long-range lists for later evaluation.
     */

    debug_gmx();

    /* Close neighbourlists */
    close_neighbor_lists(fr, bMakeQMMMnblist);

    return nns;
}

void ns_realloc_natoms(gmx_ns_t *ns, int natoms)
{
    int i;

    if (natoms > ns->nra_alloc)
    {
        ns->nra_alloc = over_alloc_dd(natoms);
        srenew(ns->bexcl, ns->nra_alloc);
        for (i = 0; i < ns->nra_alloc; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            ns->bexcl[i] = 0;
        }
    }
}

<<<<<<< HEAD
void init_ns(FILE *fplog,const t_commrec *cr,
             gmx_ns_t *ns,t_forcerec *fr,
             const gmx_mtop_t *mtop,
             matrix box)
{
    int  mt,icg,nr_in_cg,maxcg,i,j,jcg,ngid,ncg;
    t_block *cgs;
    char *ptr;
    
    /* Compute largest charge groups size (# atoms) */
    nr_in_cg=1;
    for(mt=0; mt<mtop->nmoltype; mt++) {
        cgs = &mtop->moltype[mt].cgs;
        for (icg=0; (icg < cgs->nr); icg++)
        {
            nr_in_cg=max(nr_in_cg,(int)(cgs->index[icg+1]-cgs->index[icg]));
=======
void init_ns(FILE *fplog, const t_commrec *cr,
             gmx_ns_t *ns, t_forcerec *fr,
             const gmx_mtop_t *mtop,
             matrix box)
{
    int  mt, icg, nr_in_cg, maxcg, i, j, jcg, ngid, ncg;
    t_block *cgs;
    char *ptr;

    /* Compute largest charge groups size (# atoms) */
    nr_in_cg = 1;
    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        cgs = &mtop->moltype[mt].cgs;
        for (icg = 0; (icg < cgs->nr); icg++)
        {
            nr_in_cg = max(nr_in_cg, (int)(cgs->index[icg+1]-cgs->index[icg]));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }

    /* Verify whether largest charge group is <= max cg.
<<<<<<< HEAD
     * This is determined by the type of the local exclusion type 
=======
     * This is determined by the type of the local exclusion type
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
     * Exclusions are stored in bits. (If the type is not large
     * enough, enlarge it, unsigned char -> unsigned short -> unsigned long)
     */
    maxcg = sizeof(t_excl)*8;
    if (nr_in_cg > maxcg)
    {
<<<<<<< HEAD
        gmx_fatal(FARGS,"Max #atoms in a charge group: %d > %d\n",
                  nr_in_cg,maxcg);
    }
    
    ngid = mtop->groups.grps[egcENER].nr;
    snew(ns->bExcludeAlleg,ngid);
    for(i=0; i<ngid; i++) {
        ns->bExcludeAlleg[i] = TRUE;
        for(j=0; j<ngid; j++)
=======
        gmx_fatal(FARGS, "Max #atoms in a charge group: %d > %d\n",
                  nr_in_cg, maxcg);
    }

    ngid = mtop->groups.grps[egcENER].nr;
    snew(ns->bExcludeAlleg, ngid);
    for (i = 0; i < ngid; i++)
    {
        ns->bExcludeAlleg[i] = TRUE;
        for (j = 0; j < ngid; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            if (!(fr->egp_flags[i*ngid+j] & EGP_EXCL))
            {
                ns->bExcludeAlleg[i] = FALSE;
            }
        }
    }
<<<<<<< HEAD
    
    if (fr->bGrid) {
        /* Grid search */
        ns->grid = init_grid(fplog,fr);
        init_nsgrid_lists(fr,ngid,ns);
=======

    if (fr->bGrid)
    {
        /* Grid search */
        ns->grid = init_grid(fplog, fr);
        init_nsgrid_lists(fr, ngid, ns);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        /* Simple search */
<<<<<<< HEAD
        snew(ns->ns_buf,ngid);
        for(i=0; (i<ngid); i++)
        {
            snew(ns->ns_buf[i],SHIFTS);
        }
        ncg = ncg_mtop(mtop);
        snew(ns->simple_aaj,2*ncg);
        for(jcg=0; (jcg<ncg); jcg++)
=======
        snew(ns->ns_buf, ngid);
        for (i = 0; (i < ngid); i++)
        {
            snew(ns->ns_buf[i], SHIFTS);
        }
        ncg = ncg_mtop(mtop);
        snew(ns->simple_aaj, 2*ncg);
        for (jcg = 0; (jcg < ncg); jcg++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            ns->simple_aaj[jcg]     = jcg;
            ns->simple_aaj[jcg+ncg] = jcg;
        }
    }
<<<<<<< HEAD
    
    /* Create array that determines whether or not atoms have VdW */
    snew(ns->bHaveVdW,fr->ntype);
    for(i=0; (i<fr->ntype); i++)
    {
        for(j=0; (j<fr->ntype); j++)
        {
            ns->bHaveVdW[i] = (ns->bHaveVdW[i] || 
                               (fr->bBHAM ? 
                                ((BHAMA(fr->nbfp,fr->ntype,i,j) != 0) ||
                                 (BHAMB(fr->nbfp,fr->ntype,i,j) != 0) ||
                                 (BHAMC(fr->nbfp,fr->ntype,i,j) != 0)) :
                                ((C6(fr->nbfp,fr->ntype,i,j) != 0) ||
                                 (C12(fr->nbfp,fr->ntype,i,j) != 0))));
        }
    }
    if (debug) 
        pr_bvec(debug,0,"bHaveVdW",ns->bHaveVdW,fr->ntype,TRUE);
    
    ns->nra_alloc = 0;
    ns->bexcl = NULL;
    if (!DOMAINDECOMP(cr))
    {
        /* This could be reduced with particle decomposition */
        ns_realloc_natoms(ns,mtop->natoms);
    }

    ns->nblist_initialized=FALSE;

    /* nbr list debug dump */
    {
        char *ptr=getenv("GMX_DUMP_NL");
        if (ptr)
        {
            ns->dump_nl=strtol(ptr,NULL,10);
=======

    /* Create array that determines whether or not atoms have VdW */
    snew(ns->bHaveVdW, fr->ntype);
    for (i = 0; (i < fr->ntype); i++)
    {
        for (j = 0; (j < fr->ntype); j++)
        {
            ns->bHaveVdW[i] = (ns->bHaveVdW[i] ||
                               (fr->bBHAM ?
                                ((BHAMA(fr->nbfp, fr->ntype, i, j) != 0) ||
                                 (BHAMB(fr->nbfp, fr->ntype, i, j) != 0) ||
                                 (BHAMC(fr->nbfp, fr->ntype, i, j) != 0)) :
                                ((C6(fr->nbfp, fr->ntype, i, j) != 0) ||
                                 (C12(fr->nbfp, fr->ntype, i, j) != 0))));
        }
    }
    if (debug)
    {
        pr_bvec(debug, 0, "bHaveVdW", ns->bHaveVdW, fr->ntype, TRUE);
    }

    ns->nra_alloc = 0;
    ns->bexcl     = NULL;
    if (!DOMAINDECOMP(cr))
    {
        /* This could be reduced with particle decomposition */
        ns_realloc_natoms(ns, mtop->natoms);
    }

    ns->nblist_initialized = FALSE;

    /* nbr list debug dump */
    {
        char *ptr = getenv("GMX_DUMP_NL");
        if (ptr)
        {
            ns->dump_nl = strtol(ptr, NULL, 10);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (fplog)
            {
                fprintf(fplog, "GMX_DUMP_NL = %d", ns->dump_nl);
            }
        }
        else
        {
<<<<<<< HEAD
            ns->dump_nl=0;
=======
            ns->dump_nl = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
}

<<<<<<< HEAD
			 
int search_neighbours(FILE *log,t_forcerec *fr,
                      rvec x[],matrix box,
                      gmx_localtop_t *top,
                      gmx_groups_t *groups,
                      t_commrec *cr,
                      t_nrnb *nrnb,t_mdatoms *md,
                      real lambda,real *dvdlambda,
                      gmx_grppairener_t *grppener,
                      gmx_bool bFillGrid,
                      gmx_bool bDoLongRange,
                      gmx_bool bDoForces,rvec *f)
{
    t_block  *cgs=&(top->cgs);
    rvec     box_size,grid_x0,grid_x1;
    int      i,j,m,ngid;
    real     min_size,grid_dens;
=======

int search_neighbours(FILE *log, t_forcerec *fr,
                      rvec x[], matrix box,
                      gmx_localtop_t *top,
                      gmx_groups_t *groups,
                      t_commrec *cr,
                      t_nrnb *nrnb, t_mdatoms *md,
                      real *lambda, real *dvdlambda,
                      gmx_grppairener_t *grppener,
                      gmx_bool bFillGrid,
                      gmx_bool bDoLongRangeNS,
                      gmx_bool bPadListsForKernels)
{
    t_block  *cgs = &(top->cgs);
    rvec     box_size, grid_x0, grid_x1;
    int      i, j, m, ngid;
    real     min_size, grid_dens;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    int      nsearch;
    gmx_bool     bGrid;
    char     *ptr;
    gmx_bool     *i_egp_flags;
<<<<<<< HEAD
    int      cg_start,cg_end,start,end;
=======
    int      cg_start, cg_end, start, end;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_ns_t *ns;
    t_grid   *grid;
    gmx_domdec_zones_t *dd_zones;
    put_in_list_t *put_in_list;
<<<<<<< HEAD
	
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    ns = &fr->ns;

    /* Set some local variables */
    bGrid = fr->bGrid;
<<<<<<< HEAD
    ngid = groups->grps[egcENER].nr;
    
    for(m=0; (m<DIM); m++)
    {
        box_size[m] = box[m][m];
    }
  
    if (fr->ePBC != epbcNONE)
    {
        if (sqr(fr->rlistlong) >= max_cutoff2(fr->ePBC,box))
        {
            gmx_fatal(FARGS,"One of the box vectors has become shorter than twice the cut-off length or box_yy-|box_zy| or box_zz has become smaller than the cut-off.");
        }
        if (!bGrid)
        {
            min_size = min(box_size[XX],min(box_size[YY],box_size[ZZ]));
            if (2*fr->rlistlong >= min_size)
                gmx_fatal(FARGS,"One of the box diagonal elements has become smaller than twice the cut-off length.");
        }
    }
    
    if (DOMAINDECOMP(cr))
    {
        ns_realloc_natoms(ns,cgs->index[cgs->nr]);
    }
    debug_gmx();
    
    /* Reset the neighbourlists */
    reset_neighbor_list(fr,FALSE,-1,-1);
    
    if (bGrid && bFillGrid)
    {
		
=======
    ngid  = groups->grps[egcENER].nr;

    for (m = 0; (m < DIM); m++)
    {
        box_size[m] = box[m][m];
    }

    if (fr->ePBC != epbcNONE)
    {
        if (sqr(fr->rlistlong) >= max_cutoff2(fr->ePBC, box))
        {
            gmx_fatal(FARGS, "One of the box vectors has become shorter than twice the cut-off length or box_yy-|box_zy| or box_zz has become smaller than the cut-off.");
        }
        if (!bGrid)
        {
            min_size = min(box_size[XX], min(box_size[YY], box_size[ZZ]));
            if (2*fr->rlistlong >= min_size)
            {
                gmx_fatal(FARGS, "One of the box diagonal elements has become smaller than twice the cut-off length.");
            }
        }
    }

    if (DOMAINDECOMP(cr))
    {
        ns_realloc_natoms(ns, cgs->index[cgs->nr]);
    }
    debug_gmx();

    /* Reset the neighbourlists */
    reset_neighbor_lists(fr, TRUE, TRUE);

    if (bGrid && bFillGrid)
    {

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        grid = ns->grid;
        if (DOMAINDECOMP(cr))
        {
            dd_zones = domdec_zones(cr->dd);
        }
        else
        {
            dd_zones = NULL;

<<<<<<< HEAD
            get_nsgrid_boundaries(grid,NULL,box,NULL,NULL,NULL,
                                  cgs->nr,fr->cg_cm,grid_x0,grid_x1,&grid_dens);

            grid_first(log,grid,NULL,NULL,fr->ePBC,box,grid_x0,grid_x1,
                       fr->rlistlong,grid_dens);
        }
        debug_gmx();
        
=======
            get_nsgrid_boundaries(grid->nboundeddim, box, NULL, NULL, NULL, NULL,
                                  cgs->nr, fr->cg_cm, grid_x0, grid_x1, &grid_dens);

            grid_first(log, grid, NULL, NULL, fr->ePBC, box, grid_x0, grid_x1,
                       fr->rlistlong, grid_dens);
        }
        debug_gmx();

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* Don't know why this all is... (DvdS 3/99) */
#ifndef SEGV
        start = 0;
        end   = cgs->nr;
#else
        start = fr->cg0;
        end   = (cgs->nr+1)/2;
#endif
<<<<<<< HEAD
        
        if (DOMAINDECOMP(cr))
        {
            end = cgs->nr;
            fill_grid(log,dd_zones,grid,end,-1,end,fr->cg_cm);
=======

        if (DOMAINDECOMP(cr))
        {
            end = cgs->nr;
            fill_grid(log, dd_zones, grid, end, -1, end, fr->cg_cm);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            grid->icg0 = 0;
            grid->icg1 = dd_zones->izone[dd_zones->nizone-1].cg1;
        }
        else
        {
<<<<<<< HEAD
            fill_grid(log,NULL,grid,cgs->nr,fr->cg0,fr->hcg,fr->cg_cm);
            grid->icg0 = fr->cg0;
            grid->icg1 = fr->hcg;
            debug_gmx();
            
            if (PARTDECOMP(cr))
                mv_grid(cr,grid);
            debug_gmx();
        }
        
        calc_elemnr(log,grid,start,end,cgs->nr);
        calc_ptrs(grid);
        grid_last(log,grid,start,end,cgs->nr);
        
        if (gmx_debug_at)
        {
            check_grid(debug,grid);
            print_grid(debug,grid);
=======
            fill_grid(log, NULL, grid, cgs->nr, fr->cg0, fr->hcg, fr->cg_cm);
            grid->icg0 = fr->cg0;
            grid->icg1 = fr->hcg;
            debug_gmx();

            if (PARTDECOMP(cr))
            {
                mv_grid(cr, grid);
            }
            debug_gmx();
        }

        calc_elemnr(log, grid, start, end, cgs->nr);
        calc_ptrs(grid);
        grid_last(log, grid, start, end, cgs->nr);

        if (gmx_debug_at)
        {
            check_grid(debug, grid);
            print_grid(debug, grid);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
    else if (fr->n_tpi)
    {
        /* Set the grid cell index for the test particle only.
         * The cell to cg index is not corrected, but that does not matter.
         */
<<<<<<< HEAD
        fill_grid(log,NULL,ns->grid,fr->hcg,fr->hcg-1,fr->hcg,fr->cg_cm);
    }
    debug_gmx();
    
    if (!fr->ns.bCGlist)
    {
        put_in_list = put_in_list_at;
    }
    else
    {
        put_in_list = put_in_list_cg;
=======
        fill_grid(log, NULL, ns->grid, fr->hcg, fr->hcg-1, fr->hcg, fr->cg_cm);
    }
    debug_gmx();

    if (fr->adress_type == eAdressOff)
    {
        if (!fr->ns.bCGlist)
        {
            put_in_list = put_in_list_at;
        }
        else
        {
            put_in_list = put_in_list_cg;
        }
    }
    else
    {
        put_in_list = put_in_list_adress;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    /* Do the core! */
    if (bGrid)
    {
<<<<<<< HEAD
        grid = ns->grid;
        nsearch = nsgrid_core(log,cr,fr,box,box_size,ngid,top,
                              grid,x,ns->bexcl,ns->bExcludeAlleg,
                              nrnb,md,lambda,dvdlambda,grppener,
                              put_in_list,ns->bHaveVdW,
                              bDoLongRange,bDoForces,f,
                              FALSE);
        
=======
        grid    = ns->grid;
        nsearch = nsgrid_core(log, cr, fr, box, box_size, ngid, top,
                              grid, x, ns->bexcl, ns->bExcludeAlleg,
                              nrnb, md, lambda, dvdlambda, grppener,
                              put_in_list, ns->bHaveVdW,
                              bDoLongRangeNS, FALSE);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* neighbour searching withouth QMMM! QM atoms have zero charge in
         * the classical calculation. The charge-charge interaction
         * between QM and MM atoms is handled in the QMMM core calculation
         * (see QMMM.c). The VDW however, we'd like to compute classically
         * and the QM MM atom pairs have just been put in the
         * corresponding neighbourlists. in case of QMMM we still need to
         * fill a special QMMM neighbourlist that contains all neighbours
<<<<<<< HEAD
         * of the QM atoms. If bQMMM is true, this list will now be made: 
         */
        if (fr->bQMMM && fr->qr->QMMMscheme!=eQMMMschemeoniom)
        {
            nsearch += nsgrid_core(log,cr,fr,box,box_size,ngid,top,
                                   grid,x,ns->bexcl,ns->bExcludeAlleg,
                                   nrnb,md,lambda,dvdlambda,grppener,
                                   put_in_list_qmmm,ns->bHaveVdW,
                                   bDoLongRange,bDoForces,f,
                                   TRUE);
        }
    }
    else 
    {
        nsearch = ns_simple_core(fr,top,md,box,box_size,
                                 ns->bexcl,ns->simple_aaj,
                                 ngid,ns->ns_buf,put_in_list,ns->bHaveVdW);
    }
    debug_gmx();
    
#ifdef DEBUG
    pr_nsblock(log);
#endif
    
    inc_nrnb(nrnb,eNR_NS,nsearch);
    /* inc_nrnb(nrnb,eNR_LR,fr->nlr); */
    
    return nsearch;
}

int natoms_beyond_ns_buffer(t_inputrec *ir,t_forcerec *fr,t_block *cgs,
                            matrix scale_tot,rvec *x)
{
    int  cg0,cg1,cg,a0,a1,a,i,j;
    real rint,hbuf2,scale;
    rvec *cg_cm,cgsc;
    gmx_bool bIsotropic;
    int  nBeyond;
    
    nBeyond = 0;
    
    rint = max(ir->rcoulomb,ir->rvdw);
    if (ir->rlist < rint)
    {
        gmx_fatal(FARGS,"The neighbor search buffer has negative size: %f nm",
                  ir->rlist - rint);
    }
    cg_cm = fr->cg_cm;
    
    cg0 = fr->cg0;
    cg1 = fr->hcg;
    
    if (!EI_DYNAMICS(ir->eI) || !DYNAMIC_BOX(*ir))
    {
        hbuf2 = sqr(0.5*(ir->rlist - rint));
        for(cg=cg0; cg<cg1; cg++)
        {
            a0 = cgs->index[cg];
            a1 = cgs->index[cg+1];
            for(a=a0; a<a1; a++)
            {
                if (distance2(cg_cm[cg],x[a]) > hbuf2)
=======
         * of the QM atoms. If bQMMM is true, this list will now be made:
         */
        if (fr->bQMMM && fr->qr->QMMMscheme != eQMMMschemeoniom)
        {
            nsearch += nsgrid_core(log, cr, fr, box, box_size, ngid, top,
                                   grid, x, ns->bexcl, ns->bExcludeAlleg,
                                   nrnb, md, lambda, dvdlambda, grppener,
                                   put_in_list_qmmm, ns->bHaveVdW,
                                   bDoLongRangeNS, TRUE);
        }
    }
    else
    {
        nsearch = ns_simple_core(fr, top, md, box, box_size,
                                 ns->bexcl, ns->simple_aaj,
                                 ngid, ns->ns_buf, put_in_list, ns->bHaveVdW);
    }
    debug_gmx();

#ifdef DEBUG
    pr_nsblock(log);
#endif

    inc_nrnb(nrnb, eNR_NS, nsearch);
    /* inc_nrnb(nrnb,eNR_LR,fr->nlr); */

    return nsearch;
}

int natoms_beyond_ns_buffer(t_inputrec *ir, t_forcerec *fr, t_block *cgs,
                            matrix scale_tot, rvec *x)
{
    int  cg0, cg1, cg, a0, a1, a, i, j;
    real rint, hbuf2, scale;
    rvec *cg_cm, cgsc;
    gmx_bool bIsotropic;
    int  nBeyond;

    nBeyond = 0;

    rint = max(ir->rcoulomb, ir->rvdw);
    if (ir->rlist < rint)
    {
        gmx_fatal(FARGS, "The neighbor search buffer has negative size: %f nm",
                  ir->rlist - rint);
    }
    cg_cm = fr->cg_cm;

    cg0 = fr->cg0;
    cg1 = fr->hcg;

    if (!EI_DYNAMICS(ir->eI) || !DYNAMIC_BOX(*ir))
    {
        hbuf2 = sqr(0.5*(ir->rlist - rint));
        for (cg = cg0; cg < cg1; cg++)
        {
            a0 = cgs->index[cg];
            a1 = cgs->index[cg+1];
            for (a = a0; a < a1; a++)
            {
                if (distance2(cg_cm[cg], x[a]) > hbuf2)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    nBeyond++;
                }
            }
        }
    }
    else
    {
        bIsotropic = TRUE;
<<<<<<< HEAD
        scale = scale_tot[0][0];
        for(i=1; i<DIM; i++)
=======
        scale      = scale_tot[0][0];
        for (i = 1; i < DIM; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* With anisotropic scaling, the original spherical ns volumes become
             * ellipsoids. To avoid costly transformations we use the minimum
             * eigenvalue of the scaling matrix for determining the buffer size.
             * Since the lower half is 0, the eigenvalues are the diagonal elements.
             */
<<<<<<< HEAD
            scale = min(scale,scale_tot[i][i]);
=======
            scale = min(scale, scale_tot[i][i]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (scale_tot[i][i] != scale_tot[i-1][i-1])
            {
                bIsotropic = FALSE;
            }
<<<<<<< HEAD
            for(j=0; j<i; j++)
=======
            for (j = 0; j < i; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                if (scale_tot[i][j] != 0)
                {
                    bIsotropic = FALSE;
                }
            }
        }
        hbuf2 = sqr(0.5*(scale*ir->rlist - rint));
        if (bIsotropic)
        {
<<<<<<< HEAD
            for(cg=cg0; cg<cg1; cg++)
            {
                svmul(scale,cg_cm[cg],cgsc);
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                for(a=a0; a<a1; a++)
                {
                    if (distance2(cgsc,x[a]) > hbuf2)
                    {                    
=======
            for (cg = cg0; cg < cg1; cg++)
            {
                svmul(scale, cg_cm[cg], cgsc);
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                for (a = a0; a < a1; a++)
                {
                    if (distance2(cgsc, x[a]) > hbuf2)
                    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                        nBeyond++;
                    }
                }
            }
        }
        else
        {
            /* Anistropic scaling */
<<<<<<< HEAD
            for(cg=cg0; cg<cg1; cg++)
=======
            for (cg = cg0; cg < cg1; cg++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                /* Since scale_tot contains the transpose of the scaling matrix,
                 * we need to multiply with the transpose.
                 */
<<<<<<< HEAD
                tmvmul_ur0(scale_tot,cg_cm[cg],cgsc);
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                for(a=a0; a<a1; a++)
                {
                    if (distance2(cgsc,x[a]) > hbuf2)
=======
                tmvmul_ur0(scale_tot, cg_cm[cg], cgsc);
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                for (a = a0; a < a1; a++)
                {
                    if (distance2(cgsc, x[a]) > hbuf2)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        nBeyond++;
                    }
                }
            }
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return nBeyond;
}
