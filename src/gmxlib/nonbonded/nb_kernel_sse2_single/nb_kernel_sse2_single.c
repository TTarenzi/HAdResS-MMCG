<<<<<<< HEAD
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
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
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Must come directly after config.h */
#ifdef GMX_THREAD_SHM_FDECOMP
#include <thread_mpi.h>
#endif

#include <types/simple.h>
#include <types/nrnb.h>

#include "nb_kernel_sse2_single.h"

/* Include single precision SSE intrinsics kernel headers in local directory */
#include "nb_kernel400_sse2_single.h"
#include "nb_kernel410_sse2_single.h"
#include "nb_kernel430_sse2_single.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef _MSC_VER
/* MSVC definition for __cpuid() */
#include <intrin.h>
#endif

#include "../nb_kerneltype.h"
#include "nb_kernel_sse2_single.h"

static nb_kernel_t *
kernellist_sse2_single[eNR_NBKERNEL_NR] = 
{
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel400_sse2_single,
    nb_kernel410_sse2_single,
    nb_kernel430_sse2_single
};


/* Return 0 if SSE support is present, or
 * non-zero on failure.
 */
int 
nb_kernel_sse2_single_test(FILE *                log)
{
	unsigned int level;
	unsigned int _eax,_ebx,_ecx,_edx;
	int status;
	int CPUInfo[4];
	
	if(NULL != log)
    {
		fprintf(log,"Checking CPU SSE2 support... ");
    }
        
	level = 1;
#ifdef _MSC_VER
	__cpuid(CPUInfo,1);
	
	_eax=CPUInfo[0];
	_ebx=CPUInfo[1];
	_ecx=CPUInfo[2];
	_edx=CPUInfo[3];
	
#elif defined(__x86_64__)
	/* GCC 64-bit inline asm */
	__asm__ ("push %%rbx\n\tcpuid\n\tpop %%rbx\n"                 \
			 : "=a" (_eax), "=S" (_ebx), "=c" (_ecx), "=d" (_edx) \
			 : "0" (level));
#elif defined(__i386__)
	__asm__ ("push %%ebx\n\tcpuid\n\tpop %%ebx\n"                 \
			 : "=a" (_eax), "=S" (_ebx), "=c" (_ecx), "=d" (_edx) \
			 : "0" (level));
#else
	if(NULL != log)
	{
		fprintf(log,"Don't know how to call cpuid() on this system!\n");
	}
	_eax=_ebx=_ecx=_edx=0;
#endif
        
	/* Features:                                                                                                       
	 *                                                                                                                 
	 * SSE      Bit 25 of edx should be set                                                                            
	 * SSE2     Bit 26 of edx should be set                                                                            
	 * SSE3     Bit  0 of ecx should be set                                                                            
	 * SSE4.1   Bit 19 of ecx should be set                                                                            
	 */
	status =  (_edx & (1 << 26)) != 0;

	if(NULL != log)
	{
		fprintf(log,"%s present.", (status==0) ? "not" : "");
	}
	
	/* Return SSE2 status */
	return status;
}




void
nb_kernel_setup_sse2_single(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_sse2_single_test(log) == 0)
    {
		return;
    }
	
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_sse2_single[i];
        if(p!=NULL)
		{
			list[i] = p; 
		}
    }
}    


	

	
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
 */
/*
 * Note: this file was generated by the GROMACS sse2_single kernel generator.
 */
#ifndef nb_kernel_sse2_single_h
#define nb_kernel_sse2_single_h

#include "../nb_kernel.h"

nb_kernel_t nb_kernel_ElecNone_VdwLJ_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwLJ_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwLJSh_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwLJSh_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwLJSw_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwLJSw_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecNone_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwLJ_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEw_VdwCSTab_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwLJSh_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSh_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwLJSw_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecEwSw_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwLJ_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCoul_VdwCSTab_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwLJ_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecCSTab_VdwCSTab_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecGB_VdwLJ_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecGB_VdwLJ_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecGB_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecGB_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecGB_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecGB_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSh_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwLJSw_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRFCut_VdwCSTab_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwLJ_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwNone_GeomW4W4_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomP1P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomP1P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW3P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW3P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW3W3_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW3W3_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW4P1_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW4P1_F_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW4W4_VF_sse2_single;
nb_kernel_t nb_kernel_ElecRF_VdwCSTab_GeomW4W4_F_sse2_single;


nb_kernel_info_t
    kernellist_sse2_single[] =
{
    { nb_kernel_ElecNone_VdwLJ_GeomP1P1_VF_sse2_single, "nb_kernel_ElecNone_VdwLJ_GeomP1P1_VF_sse2_single", "sse2_single", "None", "None", "LennardJones", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecNone_VdwLJ_GeomP1P1_F_sse2_single, "nb_kernel_ElecNone_VdwLJ_GeomP1P1_F_sse2_single", "sse2_single", "None", "None", "LennardJones", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecNone_VdwLJSh_GeomP1P1_VF_sse2_single, "nb_kernel_ElecNone_VdwLJSh_GeomP1P1_VF_sse2_single", "sse2_single", "None", "None", "LennardJones", "PotentialShift", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecNone_VdwLJSh_GeomP1P1_F_sse2_single, "nb_kernel_ElecNone_VdwLJSh_GeomP1P1_F_sse2_single", "sse2_single", "None", "None", "LennardJones", "PotentialShift", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecNone_VdwLJSw_GeomP1P1_VF_sse2_single, "nb_kernel_ElecNone_VdwLJSw_GeomP1P1_VF_sse2_single", "sse2_single", "None", "None", "LennardJones", "PotentialSwitch", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecNone_VdwLJSw_GeomP1P1_F_sse2_single, "nb_kernel_ElecNone_VdwLJSw_GeomP1P1_F_sse2_single", "sse2_single", "None", "None", "LennardJones", "PotentialSwitch", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecNone_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecNone_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "None", "None", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecNone_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecNone_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "None", "None", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEw_VdwLJ_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwLJ_GeomP1P1_F_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEw_VdwLJ_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwLJ_GeomW3P1_F_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEw_VdwLJ_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwLJ_GeomW3W3_F_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEw_VdwLJ_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwLJ_GeomW4P1_F_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEw_VdwLJ_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwLJ_GeomW4W4_F_sse2_single, "nb_kernel_ElecEw_VdwLJ_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "None", "LennardJones", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecEw_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEw_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEw_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEw_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEw_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecEw_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "None", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecEw_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW3P1_F_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW3W3_F_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW4P1_F_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEw_VdwCSTab_GeomW4W4_F_sse2_single, "nb_kernel_ElecEw_VdwCSTab_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "None", "CubicSplineTable", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomP1P1_F_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW3P1_F_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW3W3_F_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW4P1_F_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwLJSh_GeomW4W4_F_sse2_single, "nb_kernel_ElecEwSh_VdwLJSh_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "LennardJones", "PotentialShift", "Water4Water4", "", "Force" },
    { nb_kernel_ElecEwSh_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSh_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecEwSh_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "PotentialShift", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomP1P1_F_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW3P1_F_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW3W3_F_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW4P1_F_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwLJSw_GeomW4W4_F_sse2_single, "nb_kernel_ElecEwSw_VdwLJSw_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "LennardJones", "PotentialSwitch", "Water4Water4", "", "Force" },
    { nb_kernel_ElecEwSw_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecEwSw_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecEwSw_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "Ewald", "PotentialSwitch", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecCoul_VdwLJ_GeomP1P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomP1P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwLJ_GeomP1P1_F_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomP1P1_F_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW3P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW3P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW3P1_F_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW3P1_F_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW3W3_VF_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW3W3_VF_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW3W3_F_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW3W3_F_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW4P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW4P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW4P1_F_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW4P1_F_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW4W4_VF_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW4W4_VF_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwLJ_GeomW4W4_F_sse2_single, "nb_kernel_ElecCoul_VdwLJ_GeomW4W4_F_sse2_single", "sse2_single", "Coulomb", "None", "LennardJones", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecCoul_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecCoul_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecCoul_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecCoul_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecCoul_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecCoul_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "Coulomb", "None", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW3P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW3P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW3P1_F_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW3P1_F_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW3W3_VF_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW3W3_VF_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW3W3_F_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW3W3_F_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW4P1_VF_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW4P1_VF_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW4P1_F_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW4P1_F_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW4W4_VF_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW4W4_VF_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecCoul_VdwCSTab_GeomW4W4_F_sse2_single, "nb_kernel_ElecCoul_VdwCSTab_GeomW4W4_F_sse2_single", "sse2_single", "Coulomb", "None", "CubicSplineTable", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomP1P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomP1P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomP1P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomP1P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW3P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW3P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW3P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW3P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW3W3_VF_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW3W3_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW3W3_F_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW3W3_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW4P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW4P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW4P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW4P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW4W4_VF_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW4W4_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwLJ_GeomW4W4_F_sse2_single, "nb_kernel_ElecCSTab_VdwLJ_GeomW4W4_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "LennardJones", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecCSTab_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecCSTab_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW3P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW3P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW3P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW3P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW3W3_VF_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW3W3_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW3W3_F_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW3W3_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW4P1_VF_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW4P1_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW4P1_F_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW4P1_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW4W4_VF_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW4W4_VF_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecCSTab_VdwCSTab_GeomW4W4_F_sse2_single, "nb_kernel_ElecCSTab_VdwCSTab_GeomW4W4_F_sse2_single", "sse2_single", "CubicSplineTable", "None", "CubicSplineTable", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecGB_VdwLJ_GeomP1P1_VF_sse2_single, "nb_kernel_ElecGB_VdwLJ_GeomP1P1_VF_sse2_single", "sse2_single", "GeneralizedBorn", "None", "LennardJones", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecGB_VdwLJ_GeomP1P1_F_sse2_single, "nb_kernel_ElecGB_VdwLJ_GeomP1P1_F_sse2_single", "sse2_single", "GeneralizedBorn", "None", "LennardJones", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecGB_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecGB_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "GeneralizedBorn", "None", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecGB_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecGB_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "GeneralizedBorn", "None", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecGB_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecGB_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "GeneralizedBorn", "None", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecGB_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecGB_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "GeneralizedBorn", "None", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomP1P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW3P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW3W3_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW4P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSh_GeomW4W4_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSh_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialShift", "Water4Water4", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomP1P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW3P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW3W3_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW4P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwLJSw_GeomW4W4_F_sse2_single, "nb_kernel_ElecRFCut_VdwLJSw_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "LennardJones", "PotentialSwitch", "Water4Water4", "", "Force" },
    { nb_kernel_ElecRFCut_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecRFCut_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW3P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW3W3_F_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW4P1_F_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRFCut_VdwCSTab_GeomW4W4_F_sse2_single, "nb_kernel_ElecRFCut_VdwCSTab_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "ExactCutoff", "CubicSplineTable", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecRF_VdwLJ_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwLJ_GeomP1P1_F_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRF_VdwLJ_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwLJ_GeomW3P1_F_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRF_VdwLJ_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwLJ_GeomW3W3_F_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRF_VdwLJ_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwLJ_GeomW4P1_F_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRF_VdwLJ_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwLJ_GeomW4W4_F_sse2_single, "nb_kernel_ElecRF_VdwLJ_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "None", "LennardJones", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecRF_VdwNone_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwNone_GeomP1P1_F_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRF_VdwNone_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwNone_GeomW3P1_F_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRF_VdwNone_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwNone_GeomW3W3_F_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRF_VdwNone_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwNone_GeomW4P1_F_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRF_VdwNone_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwNone_GeomW4W4_F_sse2_single, "nb_kernel_ElecRF_VdwNone_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "None", "None", "None", "Water4Water4", "", "Force" },
    { nb_kernel_ElecRF_VdwCSTab_GeomP1P1_VF_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomP1P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "ParticleParticle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwCSTab_GeomP1P1_F_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomP1P1_F_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "ParticleParticle", "", "Force" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW3P1_VF_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW3P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water3Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW3P1_F_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW3P1_F_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water3Particle", "", "Force" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW3W3_VF_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW3W3_VF_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water3Water3", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW3W3_F_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW3W3_F_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water3Water3", "", "Force" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW4P1_VF_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW4P1_VF_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water4Particle", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW4P1_F_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW4P1_F_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water4Particle", "", "Force" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW4W4_VF_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW4W4_VF_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water4Water4", "", "PotentialAndForce" },
    { nb_kernel_ElecRF_VdwCSTab_GeomW4W4_F_sse2_single, "nb_kernel_ElecRF_VdwCSTab_GeomW4W4_F_sse2_single", "sse2_single", "ReactionField", "None", "CubicSplineTable", "None", "Water4Water4", "", "Force" }
};

int
    kernellist_sse2_single_size = sizeof(kernellist_sse2_single)/sizeof(kernellist_sse2_single[0]);

#endif
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
