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

<<<<<<< HEAD
#ifdef GMX_THREAD_SHM_FDECOMP
#include <thread_mpi.h>
#endif


=======
#ifdef GMX_THREAD_MPI
#include <thread_mpi.h>
#endif

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "txtdump.h"
#include "smalloc.h"
#include "ns.h"
#include "vec.h"
#include "maths.h"
#include "macros.h"
<<<<<<< HEAD
=======
#include "string2.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "force.h"
#include "names.h"
#include "main.h"
#include "xvgr.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "force.h"
#include "bondf.h"
#include "nrnb.h"
#include "smalloc.h"
#include "nonbonded.h"

<<<<<<< HEAD
#include "nb_kernel_c/nb_kernel_c.h"
#include "nb_kernel_adress_c/nb_kernel_c_adress.h"
=======
#include "nb_kernel.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "nb_free_energy.h"
#include "nb_generic.h"
#include "nb_generic_cg.h"
#include "nb_generic_adress.h"

<<<<<<< HEAD

/* 1,4 interactions uses kernel 330 directly */
#include "nb_kernel_c/nb_kernel330.h" 
#include "nb_kernel_adress_c/nb_kernel330_adress.h"

#ifdef GMX_PPC_ALTIVEC   
#include "nb_kernel_ppc_altivec/nb_kernel_ppc_altivec.h"
#endif

#if defined(GMX_IA32_SSE) 
#include "nb_kernel_ia32_sse/nb_kernel_ia32_sse.h"
#endif

#if defined(GMX_IA32_SSE2) 
#include "nb_kernel_ia32_sse2/nb_kernel_ia32_sse2.h"
#endif
 
#if defined(GMX_X86_64_SSE)
#include "nb_kernel_x86_64_sse/nb_kernel_x86_64_sse.h"
#endif

#if defined(GMX_X86_64_SSE2)
#include "nb_kernel_x86_64_sse2/nb_kernel_x86_64_sse2.h"
#endif

#if defined(GMX_SSE2)
#  ifdef GMX_DOUBLE
#    include "nb_kernel_sse2_double/nb_kernel_sse2_double.h"
#  else
#    include "nb_kernel_sse2_single/nb_kernel_sse2_single.h"
#  endif
#endif

#if defined(GMX_FORTRAN)
#  ifdef GMX_DOUBLE
#    include "nb_kernel_f77_double/nb_kernel_f77_double.h"
#  else
#    include "nb_kernel_f77_single/nb_kernel_f77_single.h"
#  endif
#endif

#if (defined GMX_IA64_ASM && defined GMX_DOUBLE) 
#include "nb_kernel_ia64_double/nb_kernel_ia64_double.h"
#endif

#if (defined GMX_IA64_ASM && !defined GMX_DOUBLE)
#include "nb_kernel_ia64_single/nb_kernel_ia64_single.h"
#endif

#ifdef GMX_POWER6
#include "nb_kernel_power6/nb_kernel_power6.h"
#endif

#ifdef GMX_BLUEGENE
#include "nb_kernel_bluegene/nb_kernel_bluegene.h"
#endif



enum { TABLE_NONE, TABLE_COMBINED, TABLE_COUL, TABLE_VDW, TABLE_NR };


/* Table version for each kernel.
 */
static const int 
nb_kernel_table[eNR_NBKERNEL_NR] = 
{
  TABLE_NONE,     /* kernel010 */
  TABLE_NONE,     /* kernel020 */
  TABLE_VDW,      /* kernel030 */
  TABLE_NONE,     /* kernel100 */
  TABLE_NONE,     /* kernel101 */
  TABLE_NONE,     /* kernel102 */
  TABLE_NONE,     /* kernel103 */
  TABLE_NONE,     /* kernel104 */
  TABLE_NONE,     /* kernel110 */
  TABLE_NONE,     /* kernel111 */
  TABLE_NONE,     /* kernel112 */
  TABLE_NONE,     /* kernel113 */
  TABLE_NONE,     /* kernel114 */
  TABLE_NONE,     /* kernel120 */
  TABLE_NONE,     /* kernel121 */
  TABLE_NONE,     /* kernel122 */
  TABLE_NONE,     /* kernel123 */
  TABLE_NONE,     /* kernel124 */
  TABLE_VDW,      /* kernel130 */
  TABLE_VDW,      /* kernel131 */
  TABLE_VDW,      /* kernel132 */
  TABLE_VDW,      /* kernel133 */
  TABLE_VDW,      /* kernel134 */
  TABLE_NONE,     /* kernel200 */
  TABLE_NONE,     /* kernel201 */
  TABLE_NONE,     /* kernel202 */
  TABLE_NONE,     /* kernel203 */
  TABLE_NONE,     /* kernel204 */
  TABLE_NONE,     /* kernel210 */
  TABLE_NONE,     /* kernel211 */
  TABLE_NONE,     /* kernel212 */
  TABLE_NONE,     /* kernel213 */
  TABLE_NONE,     /* kernel214 */
  TABLE_NONE,     /* kernel220 */
  TABLE_NONE,     /* kernel221 */
  TABLE_NONE,     /* kernel222 */
  TABLE_NONE,     /* kernel223 */
  TABLE_NONE,     /* kernel224 */
  TABLE_VDW,      /* kernel230 */
  TABLE_VDW,      /* kernel231 */
  TABLE_VDW,      /* kernel232 */
  TABLE_VDW,      /* kernel233 */
  TABLE_VDW,      /* kernel234 */
  TABLE_COUL,     /* kernel300 */
  TABLE_COUL,     /* kernel301 */
  TABLE_COUL,     /* kernel302 */
  TABLE_COUL,     /* kernel303 */
  TABLE_COUL,     /* kernel304 */
  TABLE_COUL,     /* kernel310 */
  TABLE_COUL,     /* kernel311 */
  TABLE_COUL,     /* kernel312 */
  TABLE_COUL,     /* kernel313 */
  TABLE_COUL,     /* kernel314 */
  TABLE_COUL,     /* kernel320 */
  TABLE_COUL,     /* kernel321 */
  TABLE_COUL,     /* kernel322 */
  TABLE_COUL,     /* kernel323 */
  TABLE_COUL,     /* kernel324 */
  TABLE_COMBINED, /* kernel330 */
  TABLE_COMBINED, /* kernel331 */
  TABLE_COMBINED, /* kernel332 */
  TABLE_COMBINED, /* kernel333 */
  TABLE_COMBINED, /* kernel334 */
  TABLE_NONE,     /* kernel400 */
  TABLE_NONE,     /* kernel410 */
  TABLE_VDW       /* kernel430 */
};



static nb_kernel_t **
nb_kernel_list = NULL;

static nb_adress_kernel_t **
nb_kernel_list_adress = NULL;

void
gmx_setup_kernels(FILE *fplog,gmx_bool bGenericKernelOnly)
{
    int i;
        
    snew(nb_kernel_list,eNR_NBKERNEL_NR);
    
    /* Note that later calls overwrite earlier, so the preferred (fastest)
     * version should be at the end. For instance, we call SSE after 3DNow.
     */
    
    for(i=0; i<eNR_NBKERNEL_NR; i++)
    {
        nb_kernel_list[i] = NULL;
    }

    // no implemented yet for h-adress
   /* if (bGenericKernelOnly)
    {*/
        return;
    /*}
	
	if(fplog)
    {
	    fprintf(fplog,"Configuring nonbonded kernels...\n");
    }
	
    nb_kernel_setup(fplog,nb_kernel_list);
    
    if(getenv("GMX_NOOPTIMIZEDKERNELS") != NULL)
    {
        return;
    }*/

    /* Setup kernels. The last called setup routine will overwrite earlier assignments,
	 * so we should e.g. test SSE3 support _after_ SSE2 support,
     * and call e.g. Fortran setup before SSE.
	 */
    
#if defined(GMX_FORTRAN) && defined(GMX_DOUBLE)   
    nb_kernel_setup_f77_double(fplog,nb_kernel_list);
#endif
	
#if defined(GMX_FORTRAN) && !defined(GMX_DOUBLE)   
    nb_kernel_setup_f77_single(fplog,nb_kernel_list);
#endif
	
#ifdef GMX_BLUEGENE
    nb_kernel_setup_bluegene(fplog,nb_kernel_list);
#endif
	
#ifdef GMX_POWER6
    nb_kernel_setup_power6(fplog,nb_kernel_list);
#endif
    
#ifdef GMX_PPC_ALTIVEC   
    nb_kernel_setup_ppc_altivec(fplog,nb_kernel_list);
#endif
	
#if defined(GMX_IA32_SSE)
    nb_kernel_setup_ia32_sse(fplog,nb_kernel_list);
#endif
	
#if defined(GMX_IA32_SSE2)
    nb_kernel_setup_ia32_sse2(fplog,nb_kernel_list);
#endif
	
#if defined(GMX_X86_64_SSE)
    nb_kernel_setup_x86_64_sse(fplog,nb_kernel_list);
#endif
	
#if defined(GMX_X86_64_SSE2)
    nb_kernel_setup_x86_64_sse2(fplog,nb_kernel_list);
#endif

#if (defined GMX_IA64_ASM && defined GMX_DOUBLE) 
    nb_kernel_setup_ia64_double(fplog,nb_kernel_list);
#endif
	
#if (defined GMX_IA64_ASM && !defined GMX_DOUBLE)
    nb_kernel_setup_ia64_single(fplog,nb_kernel_list);
#endif
	
	if(fplog)
    {
	    fprintf(fplog,"\n\n");
    }
}

void
gmx_setup_adress_kernels(FILE *fplog,gmx_bool bGenericKernelOnly) {
    int i;

    snew(nb_kernel_list_adress, eNR_NBKERNEL_NR);

    for (i = 0; i < eNR_NBKERNEL_NR; i++) {
        nb_kernel_list_adress[i] = NULL;
    }


    // Alaways using generic kernel for adress now! No implementation yet in the other kernels...
    return;
    
    /*if (bGenericKernelOnly)
    {
        return;
    }*/

    //nb_kernel_setup_adress(fplog, nb_kernel_list_adress);
}

void do_nonbonded(t_commrec *cr,t_forcerec *fr,
                  rvec x[],rvec f[],t_mdatoms *mdatoms,t_blocka *excl,
                  real egnb[],real egcoul[],real egpol[],rvec box_size,
                  t_nrnb *nrnb,real lambda,real *dvdlambda,
                  int nls,int eNL,int flags)
{
    gmx_bool            bLR,bDoForces,bForeignLambda;
	t_nblist *      nlist;
	real *          fshift;
	int             n,n0,n1,i,i0,i1,nrnb_ind,sz;
	t_nblists       *nblists;
	gmx_bool            bWater;
	nb_kernel_t *   kernelptr;
        nb_adress_kernel_t * adresskernelptr;
	FILE *          fp;
	int             fac=0;
	int             nthreads = 1;
	int             tabletype;
	int             outeriter,inneriter;
	real *          tabledata = NULL;
	gmx_gbdata_t    gbdata;
    
        gmx_bool        bCG; /* for AdresS */
        int             k;/* for AdresS */

    bLR            = (flags & GMX_DONB_LR);
    bDoForces      = (flags & GMX_DONB_FORCES);
    bForeignLambda = (flags & GMX_DONB_FOREIGNLAMBDA); 

    bCG = FALSE;  /* for AdresS */
    adresskernelptr = NULL;

	gbdata.gb_epsilon_solvent = fr->gb_epsilon_solvent;
	gbdata.epsilon_r = fr->epsilon_r;
	gbdata.gpol               = egpol;
    
    if (!fr->adress_type==eAdressOff && !bDoForces){
        gmx_fatal(FARGS,"No force kernels not implemeted for adress");
    }

    for(i=0; i < mdatoms->nalloc; i++) mdatoms->V_tot[i]=0.0;

    if(fr->bAllvsAll) 
    {
        if(fr->bGB)
        {
#if (defined GMX_SSE2 || defined GMX_X86_64_SSE || defined GMX_X86_64_SSE2 || defined GMX_IA32_SSE || defined GMX_IA32_SSE2)
# ifdef GMX_DOUBLE
            if(fr->UseOptimizedKernels)
            {
                nb_kernel_allvsallgb_sse2_double(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                                 &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else
            {
                nb_kernel_allvsallgb(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                     &outeriter,&inneriter,&fr->AllvsAll_work);        
            }
#  else /* not double */
            if(fr->UseOptimizedKernels)
            {
                nb_kernel_allvsallgb_sse2_single(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                                 &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else
            {
                nb_kernel_allvsallgb(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                     &outeriter,&inneriter,&fr->AllvsAll_work);        
            }
#  endif /* double/single alt. */
#else /* no SSE support compiled in */
            nb_kernel_allvsallgb(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                 &outeriter,&inneriter,&fr->AllvsAll_work);                    
#endif
            inc_nrnb(nrnb,eNR_NBKERNEL_ALLVSALLGB,inneriter);
        }
        else
        { 
#if (defined GMX_SSE2 || defined GMX_X86_64_SSE || defined GMX_X86_64_SSE2 || defined GMX_IA32_SSE || defined GMX_IA32_SSE2)
# ifdef GMX_DOUBLE
            if(fr->UseOptimizedKernels)
            {
                nb_kernel_allvsall_sse2_double(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                               &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else 
            {
                nb_kernel_allvsall(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                   &outeriter,&inneriter,&fr->AllvsAll_work);            
            }
            
#  else /* not double */
            if(fr->UseOptimizedKernels)
            {
                nb_kernel_allvsall_sse2_single(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                               &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else 
            {
                nb_kernel_allvsall(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                   &outeriter,&inneriter,&fr->AllvsAll_work);            
            }

#  endif /* double/single check */
#else /* No SSE2 support compiled in */
            nb_kernel_allvsall(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                               &outeriter,&inneriter,&fr->AllvsAll_work);
#endif            
            
            inc_nrnb(nrnb,eNR_NBKERNEL_ALLVSALL,inneriter);
        }
        inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,outeriter);
        return;
    }
	
    if (eNL >= 0) 
    {
		i0 = eNL;
		i1 = i0+1;
    }
    else
    {
		i0 = 0;
		i1 = eNL_NR;
	}
	
	if (nls >= 0) 
	{
		n0 = nls;
		n1 = nls+1;
	}
	else 
	{
		n0 = 0;
		n1 = fr->nnblists;
	}
	
	if(nb_kernel_list == NULL)
    {
		gmx_fatal(FARGS,"gmx_setup_kernels has not been called");
    }
  
    fshift = fr->fshift[0];
  
	for(n=n0; (n<n1); n++) 
	{
		nblists = &fr->nblists[n];
		for(i=i0; (i<i1); i++) 
		{
			outeriter = inneriter = 0;
      
			if (bLR) 
			{
				nlist = &(nblists->nlist_lr[i]);
			}
			else
			{
				nlist = &(nblists->nlist_sr[i]);
			}
			
			if (nlist->nri > 0) 
			{
				nrnb_ind = nlist->il_code;
				
				if(nrnb_ind==eNR_NBKERNEL_FREE_ENERGY)
				{
					/* generic free energy, use combined table */
					tabledata = nblists->tab.tab;
				}
				else
				{
                    if (bForeignLambda)
                    {
                        /* We don't need the non-perturbed interactions */
                        continue;
                    }

					tabletype = nb_kernel_table[nrnb_ind];
					
					/* normal kernels, not free energy */
					if (!bDoForces)
					{
						nrnb_ind += eNR_NBKERNEL_NR/2;
					}
					
					if(tabletype == TABLE_COMBINED)
					{
						tabledata = nblists->tab.tab;
					}
					else if(tabletype == TABLE_COUL)
					{
						tabledata = nblists->coultab;
					}
					else if(tabletype == TABLE_VDW)
					{
						tabledata = nblists->vdwtab;
					}
					else
					{
						tabledata = NULL;
					}
				}
				
				nlist->count = 0;

				
				if(nlist->free_energy)
				{
					if(nlist->ivdw==2)
					{
						gmx_fatal(FARGS,"Cannot do free energy Buckingham interactions.");
					}
					
					gmx_nb_free_energy_kernel(nlist->icoul,
											  nlist->ivdw,
											  nlist->nri,
											  nlist->iinr,
											  nlist->jindex,
											  nlist->jjnr,
											  nlist->shift,
											  fr->shift_vec[0],
											  fshift,
											  nlist->gid,
											  x[0],
											  f[0],
											  mdatoms->chargeA,
											  mdatoms->chargeB,
											  fr->epsfac,
											  fr->k_rf,
											  fr->c_rf,
											  fr->ewaldcoeff,
											  egcoul,
											  mdatoms->typeA,
											  mdatoms->typeB,
											  fr->ntype,
											  fr->nbfp,
											  egnb,
											  nblists->tab.scale,
											  tabledata,
											  lambda,
											  dvdlambda,
											  fr->sc_alpha,
											  fr->sc_power,
											  fr->sc_sigma6_def,
                                              fr->sc_sigma6_min,
                                              bDoForces,
											  &outeriter,
											  &inneriter);
                }
                else if (nlist->enlist == enlistCG_CG)
                {
		    if (fr->adress_type==eAdressOff){
                    /* Call the charge group based inner loop */
                       gmx_nb_generic_cg_kernel(nlist,
                                                fr,
                                                mdatoms,
                                                x[0],
                                                f[0],
                                                fshift,
                                                egcoul,
                                                egnb,
                                                nblists->tab.scale,
                                                tabledata,
                                                &outeriter,
                                                &inneriter);
		    }
		    else
		    {
                       /*gmx_nb_generic_adress_kernel(nlist,
                                                fr,
                                                mdatoms,
                                                x[0],
                                                f[0],
                                                fshift,
                                                egcoul,
                                                egnb,
                                                nblists->tab.scale,
                                                tabledata,
                                                &outeriter,
                                                &inneriter);*/
                          gmx_fatal(FARGS,"Death & horror! Adress cgcg kernel not implemented anymore.\n");

		    }
                }
                else
                {
                    /* AdresS*/
                    /* for adress we need to determine for each energy group wether it is explicit or coarse-grained */
                    if (!fr->adress_type == eAdressOff) {                        
                        bCG = FALSE;
                       /* printf("i %d n %d atom %d energygrp %d expl %d\n", i, n, nlist->iinr[0], mdatoms->cENER[nlist->iinr[0]],
                        fr->adress_group_explicit[ mdatoms->cENER[nlist->iinr[0]] ]   );
                        printf("atom 1 energygrp %d\n",mdatoms->cENER[1] );*/
                        
                        if ( !fr->adress_group_explicit[ mdatoms->cENER[nlist->iinr[0]] ] ){
                            bCG=TRUE;
                        }
                        /* If this processor has only explicit atoms (w=1)
                          skip the coarse grained force calculation. Same for
                         only coarsegrained atoms and explicit interactions.
                         Last condition is to make sure that generic kernel is not
                         skipped*/
                        if (mdatoms->pureex && bCG && nb_kernel_list[nrnb_ind] != NULL) continue;
                        if (mdatoms->purecg && !bCG && nb_kernel_list[nrnb_ind] != NULL) continue;
                        kernelptr = NULL;
                        adresskernelptr = NULL;
                    }

                    if (fr->adress_type == eAdressOff ||
                            mdatoms->pureex ||
                            mdatoms->purecg){
                        /* if we only have to calculate pure cg/ex interactions
                         we can use the faster standard gromacs kernels*/
                        kernelptr = nb_kernel_list[nrnb_ind];
                    }else{
                        /* This processor has hybrid interactions which means
                         * we have to
                         * use our own kernels. We have two kernel types: one that
                         * calculates the forces with the explicit prefactor w1*w2
                         * and one for coarse-grained with (1-w1*w2)
                         * explicit kernels are the second part of the kernel
                         *  list */
                        if (!bCG) nrnb_ind += eNR_NBKERNEL_NR/2;                      
                        adresskernelptr = nb_kernel_list_adress[nrnb_ind];
                    }
                    
                    if (kernelptr == NULL && adresskernelptr == NULL)
                     {
                        /* Call a generic nonbonded kernel */
                        
                        /* If you want to hack/test your own interactions,
                         * do it in this routine and make sure it is called
                         * by setting the environment variable GMX_NB_GENERIC.
                         */
                        if (fr->adress_type==eAdressOff){

                        gmx_nb_generic_kernel(nlist,
                                              fr,
                                              mdatoms,
                                              x[0],
                                              f[0],
                                              fshift,
                                              egcoul,
                                              egnb,
                                              nblists->tab.scale,
                                              tabledata,
                                              &outeriter,
                                              &inneriter);
                        }else /* do generic AdResS kernels (slow)*/
                        {

                            gmx_nb_generic_adress_kernel(nlist,
                                                fr,
                                                mdatoms,
                                                x[0],
                                                f[0],
                                                fshift,
                                                egcoul,
                                                egnb,
                                                nblists->tab.scale,
                                                tabledata,
                                                &outeriter,
                                                &inneriter,
                                                bCG);
                        }


                    }
                    else
                    {
                        /* Call nonbonded kernel from function pointer */
                        if (kernelptr!=NULL){
                        (*kernelptr)( &(nlist->nri),
                                      nlist->iinr,
                                      nlist->jindex,
                                      nlist->jjnr,
                                      nlist->shift,
                                      fr->shift_vec[0],
                                      fshift,
                                      nlist->gid,
                                      x[0],
                                      f[0],
                                      mdatoms->chargeA,
                                      &(fr->epsfac),
                                      &(fr->k_rf),
                                      &(fr->c_rf),
                                      egcoul,
                                      mdatoms->typeA,
                                      &(fr->ntype),
                                      fr->nbfp,
                                      egnb,
                                      &(nblists->tab.scale),
                                      tabledata,
                                      fr->invsqrta,
                                      fr->dvda,
                                      &(fr->gbtabscale),
                                      fr->gbtab.tab,
                                      &nthreads,
                                      &(nlist->count),
                                      nlist->mtx,
                                      &outeriter,
                                      &inneriter,
                                      (real *)&gbdata);
                        }else if (adresskernelptr != NULL)
                        { /* Adress kernels */
                          (*adresskernelptr)( &(nlist->nri),
                                      nlist->iinr,
                                      nlist->jindex,
                                      nlist->jjnr,
                                      nlist->shift,
                                      fr->shift_vec[0],
                                      fshift,
                                      nlist->gid,
                                      x[0],
                                      f[0],
                                      mdatoms->chargeA,
                                      &(fr->epsfac),
                                      &(fr->k_rf),
                                      &(fr->c_rf),
                                      egcoul,
                                      mdatoms->typeA,
                                      &(fr->ntype),
                                      fr->nbfp,
                                      egnb,
                                      &(nblists->tab.scale),
                                      tabledata,
                                      fr->invsqrta,
                                      fr->dvda,
                                      &(fr->gbtabscale),
                                      fr->gbtab.tab,
                                      &nthreads,
                                      &(nlist->count),
                                      nlist->mtx,
                                      &outeriter,
                                      &inneriter,
                                      fr->adress_ex_forcecap,
                                      mdatoms->wf);
                        }
                    }
                }
                
                /* Update flop accounting */
				
				/* Outer loop in kernel */
                switch (nlist->enlist) {
                case enlistATOM_ATOM:   fac =  1; break;
                case enlistSPC_ATOM:    fac =  3; break;
                case enlistSPC_SPC:     fac =  9; break;
                case enlistTIP4P_ATOM:  fac =  4; break;
                case enlistTIP4P_TIP4P: fac = 16; break;
                case enlistCG_CG:       fac =  1; break;
                }
                inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,fac*outeriter);

                /* inner loop in kernel */
                inc_nrnb(nrnb,nrnb_ind,inneriter);
=======
/* Different default (c) and accelerated interaction-specific kernels */
#include "nb_kernel_c/nb_kernel_c.h"

#if (defined GMX_CPU_ACCELERATION_X86_SSE2) && !(defined GMX_DOUBLE)
#    include "nb_kernel_sse2_single/nb_kernel_sse2_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1) && !(defined GMX_DOUBLE)
#    include "nb_kernel_sse4_1_single/nb_kernel_sse4_1_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA) && !(defined GMX_DOUBLE)
#    include "nb_kernel_avx_128_fma_single/nb_kernel_avx_128_fma_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256) && !(defined GMX_DOUBLE)
#    include "nb_kernel_avx_256_single/nb_kernel_avx_256_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE2 && defined GMX_DOUBLE)
#    include "nb_kernel_sse2_double/nb_kernel_sse2_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1 && defined GMX_DOUBLE)
#    include "nb_kernel_sse4_1_double/nb_kernel_sse4_1_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA && defined GMX_DOUBLE)
#    include "nb_kernel_avx_128_fma_double/nb_kernel_avx_128_fma_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256 && defined GMX_DOUBLE)
#    include "nb_kernel_avx_256_double/nb_kernel_avx_256_double.h"
#endif


#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t nonbonded_setup_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif
static gmx_bool            nonbonded_setup_done  = FALSE;


void
gmx_nonbonded_setup(FILE *         fplog,
                    t_forcerec *   fr,
                    gmx_bool       bGenericKernelOnly)
{
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&nonbonded_setup_mutex);
#endif
    /* Here we are guaranteed only one thread made it. */
    if (nonbonded_setup_done == FALSE)
    {
        if (bGenericKernelOnly == FALSE)
        {
            /* Add the generic kernels to the structure stored statically in nb_kernel.c */
            nb_kernel_list_add_kernels(kernellist_c, kernellist_c_size);

            if (!(fr != NULL && fr->use_cpu_acceleration == FALSE))
            {
                /* Add interaction-specific kernels for different architectures */
                /* Single precision */
#if (defined GMX_CPU_ACCELERATION_X86_SSE2) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse2_single, kernellist_sse2_single_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse4_1_single, kernellist_sse4_1_single_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_128_fma_single, kernellist_avx_128_fma_single_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_256_single, kernellist_avx_256_single_size);
#endif
                /* Double precision */
#if (defined GMX_CPU_ACCELERATION_X86_SSE2 && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse2_double, kernellist_sse2_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1 && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse4_1_double, kernellist_sse4_1_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_128_fma_double, kernellist_avx_128_fma_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256 && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_256_double, kernellist_avx_256_double_size);
#endif
                ; /* empty statement to avoid a completely empty block */
            }
        }
        /* Create a hash for faster lookups */
        nb_kernel_list_hash_init();

        nonbonded_setup_done = TRUE;
    }
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&nonbonded_setup_mutex);
#endif
}



void
gmx_nonbonded_set_kernel_pointers(FILE *log, t_nblist *nl)
{
    const char *     elec;
    const char *     elec_mod;
    const char *     vdw;
    const char *     vdw_mod;
    const char *     geom;
    const char *     other;
    const char *     vf;

    struct
    {
        const char *  arch;
        int           simd_padding_width;
    }
    arch_and_padding[] =
    {
        /* Single precision */
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256) && !(defined GMX_DOUBLE)
        { "avx_256_single", 8 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA) && !(defined GMX_DOUBLE)
        { "avx_128_fma_single", 4 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1) && !(defined GMX_DOUBLE)
        { "sse4_1_single", 4 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE2) && !(defined GMX_DOUBLE)
        { "sse2_single", 4 },
#endif
        /* Double precision */
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256 && defined GMX_DOUBLE)
        { "avx_256_double", 4 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA && defined GMX_DOUBLE)
        /* Sic. Double precision 2-way SIMD does not require neighbor list padding,
         * since the kernels execute a loop unrolled a factor 2, followed by
         * a possible single odd-element epilogue.
         */
        { "avx_128_fma_double", 1 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE2 && defined GMX_DOUBLE)
        /* No padding - see comment above */
        { "sse2_double", 1 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1 && defined GMX_DOUBLE)
        /* No padding - see comment above */
        { "sse4_1_double", 1 },
#endif
        { "c", 1 },
    };
    int              narch = asize(arch_and_padding);
    int              i;

    if (nonbonded_setup_done == FALSE)
    {
        /* We typically call this setup routine before starting timers,
         * but if that has not been done for whatever reason we do it now.
         */
        gmx_nonbonded_setup(NULL, NULL, FALSE);
    }

    /* Not used yet */
    other = "";

    nl->kernelptr_vf = NULL;
    nl->kernelptr_v  = NULL;
    nl->kernelptr_f  = NULL;

    elec     = gmx_nbkernel_elec_names[nl->ielec];
    elec_mod = eintmod_names[nl->ielecmod];
    vdw      = gmx_nbkernel_vdw_names[nl->ivdw];
    vdw_mod  = eintmod_names[nl->ivdwmod];
    geom     = gmx_nblist_geometry_names[nl->igeometry];

    if (nl->type == GMX_NBLIST_INTERACTION_ADRESS)
    {
        nl->kernelptr_vf       = (void *) gmx_nb_generic_adress_kernel;
        nl->kernelptr_f        = (void *) gmx_nb_generic_adress_kernel;
        nl->simd_padding_width = 1;
        return;
    }

    if (nl->type == GMX_NBLIST_INTERACTION_FREE_ENERGY)
    {
        nl->kernelptr_vf       = (void *) gmx_nb_free_energy_kernel;
        nl->kernelptr_f        = (void *) gmx_nb_free_energy_kernel;
        nl->simd_padding_width = 1;
    }
    else if (!gmx_strcasecmp_min(geom, "CG-CG"))
    {
        nl->kernelptr_vf       = (void *) gmx_nb_generic_cg_kernel;
        nl->kernelptr_f        = (void *) gmx_nb_generic_cg_kernel;
        nl->simd_padding_width = 1;
    }
    else
    {
        /* Try to find a specific kernel first */

        for (i = 0; i < narch && nl->kernelptr_vf == NULL; i++)
        {
            nl->kernelptr_vf       = (void *) nb_kernel_list_findkernel(log, arch_and_padding[i].arch, elec, elec_mod, vdw, vdw_mod, geom, other, "PotentialAndForce");
            nl->simd_padding_width = arch_and_padding[i].simd_padding_width;
        }
        for (i = 0; i < narch && nl->kernelptr_f == NULL; i++)
        {
            nl->kernelptr_f        = (void *) nb_kernel_list_findkernel(log, arch_and_padding[i].arch, elec, elec_mod, vdw, vdw_mod, geom, other, "Force");
            nl->simd_padding_width = arch_and_padding[i].simd_padding_width;

            /* If there is not force-only optimized kernel, is there a potential & force one? */
            if (nl->kernelptr_f == NULL)
            {
                nl->kernelptr_f        = (void *) nb_kernel_list_findkernel(NULL, arch_and_padding[i].arch, elec, elec_mod, vdw, vdw_mod, geom, other, "PotentialAndForce");
                nl->simd_padding_width = arch_and_padding[i].simd_padding_width;
            }
        }

        /* Give up, pick a generic one instead */
        if (nl->kernelptr_vf == NULL)
        {
            nl->kernelptr_vf       = (void *) gmx_nb_generic_kernel;
            nl->kernelptr_f        = (void *) gmx_nb_generic_kernel;
            nl->simd_padding_width = 1;
            if (debug)
            {
                fprintf(debug,
                        "WARNING - Slow generic NB kernel used for neighborlist with\n"
                        "    Elec: '%s', Modifier: '%s'\n"
                        "    Vdw:  '%s', Modifier: '%s'\n"
                        "    Geom: '%s', Other: '%s'\n\n",
                        elec, elec_mod, vdw, vdw_mod, geom, other);
            }
        }
    }

    return;
}

void do_nonbonded(t_commrec *cr, t_forcerec *fr,
                  rvec x[], rvec f_shortrange[], rvec f_longrange[], t_mdatoms *mdatoms, t_blocka *excl,
                  gmx_grppairener_t *grppener, rvec box_size,
                  t_nrnb *nrnb, real *lambda, real *dvdl,
                  int nls, int eNL, int flags)
{
    t_nblist *        nlist;
    int               n, n0, n1, i, i0, i1, sz, range;
    t_nblists *       nblists;
    nb_kernel_data_t  kernel_data;
    nb_kernel_t *     kernelptr = NULL;
    rvec *            f;

    kernel_data.flags                   = flags;
    kernel_data.exclusions              = excl;
    kernel_data.lambda                  = lambda;
    kernel_data.dvdl                    = dvdl;

    if (fr->bAllvsAll)
    {
        return;
    }

    if (eNL >= 0)
    {
        i0 = eNL;
        i1 = i0+1;
    }
    else
    {
        i0 = 0;
        i1 = eNL_NR;
    }

    if (nls >= 0)
    {
        n0 = nls;
        n1 = nls+1;
    }
    else
    {
        n0 = 0;
        n1 = fr->nnblists;
    }

    for (n = n0; (n < n1); n++)
    {
        nblists = &fr->nblists[n];

        kernel_data.table_elec              = &nblists->table_elec;
        kernel_data.table_vdw               = &nblists->table_vdw;
        kernel_data.table_elec_vdw          = &nblists->table_elec_vdw;

        for (range = 0; range < 2; range++)
        {
            /* Are we doing short/long-range? */
            if (range == 0)
            {
                /* Short-range */
                if (!(flags & GMX_NONBONDED_DO_SR))
                {
                    continue;
                }
                kernel_data.energygrp_elec          = grppener->ener[egCOULSR];
                kernel_data.energygrp_vdw           = grppener->ener[fr->bBHAM ? egBHAMSR : egLJSR];
                kernel_data.energygrp_polarization  = grppener->ener[egGB];
                nlist = nblists->nlist_sr;
                f                                   = f_shortrange;
            }
            else if (range == 1)
            {
                /* Long-range */
                if (!(flags & GMX_NONBONDED_DO_LR))
                {
                    continue;
                }
                kernel_data.energygrp_elec          = grppener->ener[egCOULLR];
                kernel_data.energygrp_vdw           = grppener->ener[fr->bBHAM ? egBHAMLR : egLJLR];
                kernel_data.energygrp_polarization  = grppener->ener[egGB];
                nlist = nblists->nlist_lr;
                f                                   = f_longrange;
            }

            for (i = i0; (i < i1); i++)
            {
                if (nlist[i].nri > 0)
                {
                    if (flags & GMX_NONBONDED_DO_POTENTIAL)
                    {
                        /* Potential and force */
                        kernelptr = (nb_kernel_t *)nlist[i].kernelptr_vf;
                    }
                    else
                    {
                        /* Force only, no potential */
                        kernelptr = (nb_kernel_t *)nlist[i].kernelptr_f;
                    }

                    if (nlist[i].type != GMX_NBLIST_INTERACTION_FREE_ENERGY && (flags & GMX_NONBONDED_DO_FOREIGNLAMBDA))
                    {
                        /* We don't need the non-perturbed interactions */
                        continue;
                    }
                    (*kernelptr)(&(nlist[i]), x, f, fr, mdatoms, &kernel_data, nrnb);
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
}

<<<<<<< HEAD

real 
do_listed_vdw_q(int ftype,int nbonds,
                const t_iatom iatoms[],const t_iparams iparams[],
                const rvec x[],rvec f[],rvec fshift[],
                const t_pbc *pbc,const t_graph *g,
                real lambda,real *dvdlambda,
                const t_mdatoms *md,
                const t_forcerec *fr,gmx_grppairener_t *grppener,
                int *global_atom_index)
{
    static    gmx_bool bWarn=FALSE;
    real      eps,r2,*tab,rtab2=0;
    rvec      dx,x14[2],f14[2];
    int       i,ai,aj,itype;
    int       typeA[2]={0,0},typeB[2]={0,1};
    real      chargeA[2]={0,0},chargeB[2];
    int       gid,shift_vir,shift_f;
    int       j_index[] = { 0, 1 };
    int       i0=0,i1=1,i2=2;
    ivec      dt;
    int       outeriter,inneriter;
    int       nthreads = 1;
    int       count;
    real      krf,crf,tabscale;
    int       ntype=0;
    real      *nbfp=NULL;
    real      *egnb=NULL,*egcoul=NULL;
    t_nblist  tmplist;
    int       icoul,ivdw;
    gmx_bool      bMolPBC,bFreeEnergy;
    
    gmx_bool      bCG; /* AdResS*/
    real      wf14[2]={0,0}; /* AdResS*/
   
#if GMX_THREAD_SHM_FDECOMP
    pthread_mutex_t mtx;
#else
    void *    mtx = NULL;
#endif

    
#if GMX_THREAD_SHM_FDECOMP
    pthread_mutex_initialize(&mtx);
#endif

    bMolPBC = fr->bMolPBC;

    switch (ftype) {
    case F_LJ14:
    case F_LJC14_Q:
        eps = fr->epsfac*fr->fudgeQQ;
        ntype  = 1;
        egnb   = grppener->ener[egLJ14];
        egcoul = grppener->ener[egCOUL14];
        break;
    case F_LJC_PAIRS_NB:
        eps = fr->epsfac;
        ntype  = 1;
        egnb   = grppener->ener[egLJSR];
        egcoul = grppener->ener[egCOULSR];
        break;
    default:
        gmx_fatal(FARGS,"Unknown function type %d in do_nonbonded14",
                  ftype);
    }
    tab = fr->tab14.tab;
    rtab2 = sqr(fr->tab14.r);
    tabscale = fr->tab14.scale;

    krf = fr->k_rf;
    crf = fr->c_rf;

    /* Determine the values for icoul/ivdw. */
    if (fr->bEwald) {
        icoul = 1;
    } 
    else if(fr->bcoultab)
    {
        icoul = 3;
    }
    else if(fr->eeltype == eelRF_NEC)
    {
        icoul = 2;
    }
    else 
    {
        icoul = 1;
    }
    
    if(fr->bvdwtab)
    {
        ivdw = 3;
    }
    else if(fr->bBHAM)
    {
        ivdw = 2;
    }
    else 
    {
        ivdw = 1;
    }
    
    
    bCG = FALSE; /*Adres*/
    /* We don't do SSE or altivec here, due to large overhead for 4-fold 
     * unrolling on short lists 
     */
    
    bFreeEnergy = FALSE;
    for(i=0; (i<nbonds); ) 
=======
static void
nb_listed_warning_rlimit(const rvec *x, int ai, int aj, int * global_atom_index, real r, real rlimit)
{
    gmx_warning("Listed nonbonded interaction between particles %d and %d\n"
                "at distance %.3f which is larger than the table limit %.3f nm.\n\n"
                "This is likely either a 1,4 interaction, or a listed interaction inside\n"
                "a smaller molecule you are decoupling during a free energy calculation.\n"
                "Since interactions at distances beyond the table cannot be computed,\n"
                "they are skipped until they are inside the table limit again. You will\n"
                "only see this message once, even if it occurs for several interactions.\n\n"
                "IMPORTANT: This should not happen in a stable simulation, so there is\n"
                "probably something wrong with your system. Only change the table-extension\n"
                "distance in the mdp file if you are really sure that is the reason.\n",
                glatnr(global_atom_index, ai), glatnr(global_atom_index, aj), r, rlimit);

    if (debug)
    {
        fprintf(debug,
                "%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. Ignored\n",
                x[ai][XX], x[ai][YY], x[ai][ZZ], x[aj][XX], x[aj][YY], x[aj][ZZ],
                glatnr(global_atom_index, ai), glatnr(global_atom_index, aj), r);
    }
}



/* This might logically belong better in the nb_generic.c module, but it is only
 * used in do_nonbonded_listed(), and we want it to be inlined there to avoid an
 * extra functional call for every single pair listed in the topology.
 */
static real
nb_evaluate_single(real r2, real tabscale, real *vftab,
                   real qq, real c6, real c12, real *velec, real *vvdw)
{
    real       rinv, r, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VVe, FFe, VVd, FFd, VVr, FFr, fscal;
    int        ntab;

    /* Do the tabulated interactions - first table lookup */
    rinv             = gmx_invsqrt(r2);
    r                = r2*rinv;
    rtab             = r*tabscale;
    ntab             = rtab;
    eps              = rtab-ntab;
    eps2             = eps*eps;
    ntab             = 12*ntab;
    /* Electrostatics */
    Y                = vftab[ntab];
    F                = vftab[ntab+1];
    Geps             = eps*vftab[ntab+2];
    Heps2            = eps2*vftab[ntab+3];
    Fp               = F+Geps+Heps2;
    VVe              = Y+eps*Fp;
    FFe              = Fp+Geps+2.0*Heps2;
    /* Dispersion */
    Y                = vftab[ntab+4];
    F                = vftab[ntab+5];
    Geps             = eps*vftab[ntab+6];
    Heps2            = eps2*vftab[ntab+7];
    Fp               = F+Geps+Heps2;
    VVd              = Y+eps*Fp;
    FFd              = Fp+Geps+2.0*Heps2;
    /* Repulsion */
    Y                = vftab[ntab+8];
    F                = vftab[ntab+9];
    Geps             = eps*vftab[ntab+10];
    Heps2            = eps2*vftab[ntab+11];
    Fp               = F+Geps+Heps2;
    VVr              = Y+eps*Fp;
    FFr              = Fp+Geps+2.0*Heps2;

    *velec           = qq*VVe;
    *vvdw            = c6*VVd+c12*VVr;

    fscal            = -(qq*FFe+c6*FFd+c12*FFr)*tabscale*rinv;

    return fscal;
}


real
do_nonbonded_listed(int ftype, int nbonds,
                    const t_iatom iatoms[], const t_iparams iparams[],
                    const rvec x[], rvec f[], rvec fshift[],
                    const t_pbc *pbc, const t_graph *g,
                    real *lambda, real *dvdl,
                    const t_mdatoms *md,
                    const t_forcerec *fr, gmx_grppairener_t *grppener,
                    int *global_atom_index)
{
    int              ielec, ivdw;
    real             qq, c6, c12;
    rvec             dx;
    ivec             dt;
    int              i, j, itype, ai, aj, gid;
    int              fshift_index;
    real             r2, rinv;
    real             fscal, velec, vvdw;
    real *           energygrp_elec;
    real *           energygrp_vdw;
    static gmx_bool  warned_rlimit = FALSE;
    /* Free energy stuff */
    gmx_bool         bFreeEnergy;
    real             LFC[2], LFV[2], DLF[2], lfac_coul[2], lfac_vdw[2], dlfac_coul[2], dlfac_vdw[2];
    real             qqB, c6B, c12B, sigma2_def, sigma2_min;


    switch (ftype)
    {
        case F_LJ14:
        case F_LJC14_Q:
            energygrp_elec = grppener->ener[egCOUL14];
            energygrp_vdw  = grppener->ener[egLJ14];
            break;
        case F_LJC_PAIRS_NB:
            energygrp_elec = grppener->ener[egCOULSR];
            energygrp_vdw  = grppener->ener[egLJSR];
            break;
        default:
            energygrp_elec = NULL; /* Keep compiler happy */
            energygrp_vdw  = NULL; /* Keep compiler happy */
            gmx_fatal(FARGS, "Unknown function type %d in do_nonbonded14", ftype);
            break;
    }

    if (fr->efep != efepNO)
    {
        /* Lambda factor for state A=1-lambda and B=lambda */
        LFC[0] = 1.0 - lambda[efptCOUL];
        LFV[0] = 1.0 - lambda[efptVDW];
        LFC[1] = lambda[efptCOUL];
        LFV[1] = lambda[efptVDW];

        /*derivative of the lambda factor for state A and B */
        DLF[0] = -1;
        DLF[1] = 1;

        /* precalculate */
        sigma2_def = pow(fr->sc_sigma6_def, 1.0/3.0);
        sigma2_min = pow(fr->sc_sigma6_min, 1.0/3.0);

        for (i = 0; i < 2; i++)
        {
            lfac_coul[i]  = (fr->sc_power == 2 ? (1-LFC[i])*(1-LFC[i]) : (1-LFC[i]));
            dlfac_coul[i] = DLF[i]*fr->sc_power/fr->sc_r_power*(fr->sc_power == 2 ? (1-LFC[i]) : 1);
            lfac_vdw[i]   = (fr->sc_power == 2 ? (1-LFV[i])*(1-LFV[i]) : (1-LFV[i]));
            dlfac_vdw[i]  = DLF[i]*fr->sc_power/fr->sc_r_power*(fr->sc_power == 2 ? (1-LFV[i]) : 1);
        }
    }
    else
    {
        sigma2_min = sigma2_def = 0;
    }

    bFreeEnergy = FALSE;
    for (i = 0; (i < nbonds); )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        itype = iatoms[i++];
        ai    = iatoms[i++];
        aj    = iatoms[i++];
<<<<<<< HEAD
        gid   = GID(md->cENER[ai],md->cENER[aj],md->nenergrp);
        
        if (!fr->adress_type == eAdressOff) {
            if (fr->adress_group_explicit[md->cENER[ai]] != fr->adress_group_explicit[md->cENER[aj]]){
                /*exclude cg-ex interaction*/
                continue;
            }           
            bCG = !fr->adress_group_explicit[md->cENER[ai]];
            wf14[0] = md->wf[ai];
            wf14[1] = md->wf[aj];
        }
        switch (ftype) {
        case F_LJ14:
            bFreeEnergy =
                (fr->efep != efepNO &&
                 ((md->nPerturbed && (md->bPerturbed[ai] || md->bPerturbed[aj])) ||
                  iparams[itype].lj14.c6A != iparams[itype].lj14.c6B ||
                  iparams[itype].lj14.c12A != iparams[itype].lj14.c12B));
            chargeA[0] = md->chargeA[ai];
            chargeA[1] = md->chargeA[aj];
            nbfp = (real *)&(iparams[itype].lj14.c6A);
            break;
        case F_LJC14_Q:
            eps = fr->epsfac*iparams[itype].ljc14.fqq;
            chargeA[0] = iparams[itype].ljc14.qi;
            chargeA[1] = iparams[itype].ljc14.qj;
            nbfp = (real *)&(iparams[itype].ljc14.c6);
            break;
        case F_LJC_PAIRS_NB:
            chargeA[0] = iparams[itype].ljcnb.qi;
            chargeA[1] = iparams[itype].ljcnb.qj;
            nbfp = (real *)&(iparams[itype].ljcnb.c6);
            break;
        }
        
        if (!bMolPBC) 
        {
            /* This is a bonded interaction, atoms are in the same box */
            shift_f = CENTRAL;
            r2 = distance2(x[ai],x[aj]);
        }
        else 
        {
            /* Apply full periodic boundary conditions */
            shift_f = pbc_dx_aiuc(pbc,x[ai],x[aj],dx);
            r2 = norm2(dx);
        }

        if (r2 >= rtab2) 
        {
            if (!bWarn) 
            {
                fprintf(stderr,"Warning: 1-4 interaction between %d and %d "
                        "at distance %.3f which is larger than the 1-4 table size %.3f nm\n", 
			glatnr(global_atom_index,ai),
			glatnr(global_atom_index,aj),
			sqrt(r2), sqrt(rtab2));
                fprintf(stderr,"These are ignored for the rest of the simulation\n");
                fprintf(stderr,"This usually means your system is exploding,\n"
                        "if not, you should increase table-extension in your mdp file\n"
                        "or with user tables increase the table size\n");
                bWarn = TRUE;
            }
            if (debug) 
	      fprintf(debug,"%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. Ignored\n",
		      x[ai][XX],x[ai][YY],x[ai][ZZ],
		      x[aj][XX],x[aj][YY],x[aj][ZZ],
		      glatnr(global_atom_index,ai),
		      glatnr(global_atom_index,aj),
		      sqrt(r2));
        }
        else 
        {
            copy_rvec(x[ai],x14[0]);
            copy_rvec(x[aj],x14[1]);
            clear_rvec(f14[0]);
            clear_rvec(f14[1]);
#ifdef DEBUG
            fprintf(debug,"LJ14: grp-i=%2d, grp-j=%2d, ngrp=%2d, GID=%d\n",
                    md->cENER[ai],md->cENER[aj],md->nenergrp,gid);
#endif
            
	    outeriter = inneriter = count = 0;
	    if (bFreeEnergy)
        {
            chargeB[0] = md->chargeB[ai];
            chargeB[1] = md->chargeB[aj];
            /* We pass &(iparams[itype].lj14.c6A) as LJ parameter matrix
             * to the innerloops.
             * Here we use that the LJ-14 parameters are stored in iparams
             * as c6A,c12A,c6B,c12B, which are referenced correctly
             * in the innerloops if we assign type combinations 0-0 and 0-1
             * to atom pair ai-aj in topologies A and B respectively.
             */
            if(ivdw==2)
            {
                gmx_fatal(FARGS,"Cannot do free energy Buckingham interactions.");
            }
            count = 0;
            gmx_nb_free_energy_kernel(icoul,
                                      ivdw,
                                      i1,
                                      &i0,
                                      j_index,
                                      &i1,
                                      &shift_f,
                                      fr->shift_vec[0],
                                      fshift[0],
                                      &gid,
                                      x14[0],
                                      f14[0],
                                      chargeA,
                                      chargeB,
                                      eps,
                                      krf,
                                      crf,
                                      fr->ewaldcoeff,
                                      egcoul,
                                      typeA,
                                      typeB,
                                      ntype,
                                      nbfp,
                                      egnb,
                                      tabscale,
                                      tab,
                                      lambda,
                                      dvdlambda,
                                      fr->sc_alpha,
                                      fr->sc_power,
                                      fr->sc_sigma6_def,
                                      fr->sc_sigma6_min,
                                      TRUE,
                                      &outeriter,
                                      &inneriter);
        }
        else 
        { 
          if (fr->adress_type==eAdressOff || !fr->adress_do_hybridpairs){
            /* Not perturbed - call kernel 330 */
            nb_kernel330
                ( &i1,
                  &i0,
                  j_index,
                  &i1,
                  &shift_f,
                  fr->shift_vec[0],
                  fshift[0],
                  &gid,
                  x14[0],
                  f14[0],
                  chargeA,
                  &eps,
                  &krf,
                  &crf,
                  egcoul,
                  typeA,
                  &ntype,
                  nbfp,
                  egnb,
                  &tabscale,
                  tab,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  &nthreads,
                  &count,
                  (void *)&mtx,
                  &outeriter,
                  &inneriter,
                  NULL);                
                } else {
                    if (bCG) {
                        nb_kernel330_adress_cg(&i1,
                                &i0,
                                j_index,
                                &i1,
                                &shift_f,
                                fr->shift_vec[0],
                                fshift[0],
                                &gid,
                                x14[0],
                                f14[0],
                                chargeA,
                                &eps,
                                &krf,
                                &crf,
                                egcoul,
                                typeA,
                                &ntype,
                                nbfp,
                                egnb,
                                &tabscale,
                                tab,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                &nthreads,
                                &count,
                                (void *) &mtx,
                                &outeriter,
                                &inneriter,
                                fr->adress_ex_forcecap,
                                wf14);
                    } else {
                        nb_kernel330_adress_ex(&i1,
                                &i0,
                                j_index,
                                &i1,
                                &shift_f,
                                fr->shift_vec[0],
                                fshift[0],
                                &gid,
                                x14[0],
                                f14[0],
                                chargeA,
                                &eps,
                                &krf,
                                &crf,
                                egcoul,
                                typeA,
                                &ntype,
                                nbfp,
                                egnb,
                                &tabscale,
                                tab,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                &nthreads,
                                &count,
                                (void *) &mtx,
                                &outeriter,
                                &inneriter,
                                fr->adress_ex_forcecap,
                                wf14);
                    }

                }
            }
        
        /* Add the forces */
        rvec_inc(f[ai],f14[0]);
        rvec_dec(f[aj],f14[0]);
        
        if (g) 
        {
            /* Correct the shift forces using the graph */
            ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);    
            shift_vir = IVEC2IS(dt);
            rvec_inc(fshift[shift_vir],f14[0]);
            rvec_dec(fshift[CENTRAL],f14[0]);
        }
        
	    /* flops: eNR_KERNEL_OUTER + eNR_KERNEL330 + 12 */
        }
    }
    return 0.0;
}


void
adress_drift_term(FILE *               fplog,
                          int                  cg0,
                          int                  cg1,
                          int                  cg1home,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc, rvec f[], real *engdelta)
{
    int            icg,k,k0,k1,d, i;
    real           nrcg,inv_ncg,mtot,inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr,adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;
    real *         wfprime;
    rvec dx;
    real r;
    real V, fscal;
    rvec fdrift;
    real Mtot;
    real massfac;
    int binlambda;
    real *dhdl, *cgdens;
    real force_cap;
    snew(dhdl, fr->adress_dhdlbins);
    snew(cgdens, fr->adress_dhdlbins);
    
    ref = &(fr->adress_refs);
    
    
    *engdelta=0.0;
    clear_rvec(fdrift);
    cgindex = cgs->index;
    force_cap = fr->adress_ex_forcecap;

    for (i = 0; i < fr->adress_dhdlbins; i++) {
         dhdl[i]=0.0;
         cgdens[i]=0.0;
     }

    if (fr->adress_onthefly_TI) {
        for (icg = cg0; (icg < cg1home); icg++) {
            k0 = cgindex[icg];
            k1 = cgindex[icg + 1];
            binlambda=(int)floor(((real)mdatoms->wf[k0]*(real)fr->adress_dhdlbins));
            if (!fr->adress_group_explicit[mdatoms->cENER[k0]] && fr->adress_cgdens[binlambda] > 0) {
                //f[k0][0] += fr->adress_fcorr[binlambda]/fr->adress_fcorr_count*mdatoms->wfprime[k0];
                f[k0][0] += fr->adress_dhdl[binlambda]/fr->adress_cgdens[binlambda]*mdatoms->wfprime[k0];
                // todo more general geometries...
            }

        }
    }


    for (icg = cg0; (icg < cg1); icg++) {
        k0 = cgindex[icg];
        k1 = cgindex[icg + 1];
        
        if (pbc)
        {
            pbc_dx(pbc,(*ref),x[k0],dx);
        }
        else
        {
            rvec_sub((*ref),x[k0],dx);
        }


        Mtot=0;
        fscal=0;
        V=0;
        for (k = k0; (k < k1); k++) {
            V+=0.5*mdatoms->V_tot[k];
            Mtot+=mdatoms->massT[k];
        }
        if (!fr->adress_group_explicit[mdatoms->cENER[k0]]){
            /* has to be plus because F=-(Vat-Vcg-DeltaV) and we will flip sign of fscal below if coarse-grained*/    
            if(icg<cg1home){
                V+=fr->adress_deltaU;
                  /* only count energy for molecules living on this node*/
                *engdelta-=fr->adress_deltaU*mdatoms->wf[k0];
             }
        }
        fscal=V*mdatoms->wfprime[k0];
        if (fr->adress_group_explicit[mdatoms->cENER[k0]]){
            fscal*=-1.0;
        }
        if (mdatoms->wf[k0]> 0.0 && mdatoms->wf[k0] < 1.0 && icg<cg1home){
            binlambda=(int)floor(((real)mdatoms->wf[k0]*(real)fr->adress_dhdlbins));
            
            if (!fr->adress_group_explicit[mdatoms->cENER[k0]]){
                cgdens[binlambda]+=1;
                dhdl[binlambda]+=V;  // + Vcg
            }else{
                dhdl[binlambda]-=V; // - Vat
            }
        }

        if (force_cap > 0 && (fabs(fscal) > force_cap)) {
            fscal = force_cap * fscal / fabs(fscal);
        }
        if (fr->adress_type == eAdressXSplit ){
            fdrift[0]=fscal;
            fdrift[1]=0.0;
            fdrift[2]=0.0;
        }else if(fr->adress_type == eAdressSphere){
            //gmx_fatal(FARGS,"Deine mutti hat spherical case not implemented!");
            r=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
            fdrift[0]=fscal*dx[0]/r;
            fdrift[1]=fscal*dx[1]/r;
            fdrift[2]=fscal*dx[2]/r;
        }

        massfac=1.0;


        for (k = k0; (k < k1); k++) {
            if (fr->adress_group_explicit[mdatoms->cENER[k]]) {
                massfac = mdatoms->massT[k] / Mtot;
              
            }

              /*if (fr->adress_group_explicit[mdatoms->cENER[k]]) {
                    printf("AA k0 %d V %g fdriftx %g massfac %g cener %d Vi %g\n", k0, V, fdrift[0], massfac, mdatoms->cENER[k], mdatoms->V_tot[k]);
                } else {
                    printf("CG k0 %d V %g fdriftx %g massfac %g cener %d Vi %g\n", k0, V, fdrift[0], massfac, mdatoms->cENER[k], mdatoms->V_tot[k]);
                }*/
        
            /* Attention: this has to be -Force because the forces in the kernel are actually F=grad V*/
            f[k][0] -= fdrift[0] * massfac;
            f[k][1] -= fdrift[1] * massfac;
            f[k][2] -= fdrift[2] * massfac;
        }

    }



     for (i = 0; i < fr->adress_dhdlbins; i++) {
           if (fr->adress_cgdens[i] > 0){
                fr->adress_fcorr[i]+=dhdl[i]/cgdens[i];
           }else{
               fr->adress_fcorr[i]+=0.0;
           }
           fr->adress_dhdl[i]+=dhdl[i];
           fr->adress_cgdens[i]+=cgdens[i];
     }
    fr->adress_fcorr_count++;

    sfree(dhdl);
    sfree(cgdens);

}

=======
        gid   = GID(md->cENER[ai], md->cENER[aj], md->nenergrp);

        /* Get parameters */
        switch (ftype)
        {
            case F_LJ14:
                bFreeEnergy =
                    (fr->efep != efepNO &&
                     ((md->nPerturbed && (md->bPerturbed[ai] || md->bPerturbed[aj])) ||
                      iparams[itype].lj14.c6A != iparams[itype].lj14.c6B ||
                      iparams[itype].lj14.c12A != iparams[itype].lj14.c12B));
                qq               = md->chargeA[ai]*md->chargeA[aj]*fr->epsfac*fr->fudgeQQ;
                c6               = iparams[itype].lj14.c6A;
                c12              = iparams[itype].lj14.c12A;
                break;
            case F_LJC14_Q:
                qq               = iparams[itype].ljc14.qi*iparams[itype].ljc14.qj*fr->epsfac*iparams[itype].ljc14.fqq;
                c6               = iparams[itype].ljc14.c6;
                c12              = iparams[itype].ljc14.c12;
                break;
            case F_LJC_PAIRS_NB:
                qq               = iparams[itype].ljcnb.qi*iparams[itype].ljcnb.qj*fr->epsfac;
                c6               = iparams[itype].ljcnb.c6;
                c12              = iparams[itype].ljcnb.c12;
                break;
            default:
                /* Cannot happen since we called gmx_fatal() above in this case */
                qq = c6 = c12 = 0; /* Keep compiler happy */
                break;
        }

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0;
        c12 *= 12.0;

        /* Do we need to apply full periodic boundary conditions? */
        if (fr->bMolPBC == TRUE)
        {
            fshift_index = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            fshift_index = CENTRAL;
            rvec_sub(x[ai], x[aj], dx);
        }
        r2           = norm2(dx);

        if (r2 >= fr->tab14.r*fr->tab14.r)
        {
            if (warned_rlimit == FALSE)
            {
                nb_listed_warning_rlimit(x, ai, aj, global_atom_index, sqrt(r2), fr->tab14.r);
                warned_rlimit = TRUE;
            }
            continue;
        }

        if (bFreeEnergy)
        {
            /* Currently free energy is only supported for F_LJ14, so no need to check for that if we got here */
            qqB              = md->chargeB[ai]*md->chargeB[aj]*fr->epsfac*fr->fudgeQQ;
            c6B              = iparams[itype].lj14.c6B*6.0;
            c12B             = iparams[itype].lj14.c12B*12.0;

            fscal            = nb_free_energy_evaluate_single(r2, fr->sc_r_power, fr->sc_alphacoul, fr->sc_alphavdw,
                                                              fr->tab14.scale, fr->tab14.data, qq, c6, c12, qqB, c6B, c12B,
                                                              LFC, LFV, DLF, lfac_coul, lfac_vdw, dlfac_coul, dlfac_vdw,
                                                              fr->sc_sigma6_def, fr->sc_sigma6_min, sigma2_def, sigma2_min, &velec, &vvdw, dvdl);
        }
        else
        {
            /* Evaluate tabulated interaction without free energy */
            fscal            = nb_evaluate_single(r2, fr->tab14.scale, fr->tab14.data, qq, c6, c12, &velec, &vvdw);
        }

        energygrp_elec[gid]  += velec;
        energygrp_vdw[gid]   += vvdw;
        svmul(fscal, dx, dx);

        /* Add the forces */
        rvec_inc(f[ai], dx);
        rvec_dec(f[aj], dx);

        if (g)
        {
            /* Correct the shift forces using the graph */
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            fshift_index = IVEC2IS(dt);
        }
        if (fshift_index != CENTRAL)
        {
            rvec_inc(fshift[fshift_index], dx);
            rvec_dec(fshift[CENTRAL], dx);
        }
    }
    return 0.0;
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
