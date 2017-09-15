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
#ifndef _NB_KERNEL_SSE2_SINGLE_H_
#define _NB_KERNEL_SSE2_SINGLE_H_

/*! \file  nb_kernel_sse2_single.h
 *  \brief SSE2-intrinsics optimized level2 nonbonded kernels.
 *
 *  \internal
 */

#include <stdio.h>

#include <types/simple.h>

#include "../nb_kerneltype.h"
#include "nb_kernel_allvsall_sse2_single.h"
#include "nb_kernel_allvsallgb_sse2_single.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

void
nb_kernel_setup_sse2_single(FILE *log,nb_kernel_t **list);

#ifdef __cplusplus
}
#endif

#endif /* _NB_KERNEL_SSE2_SINGLE_H_ */
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2001-2012, The GROMACS Development Team
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
#ifndef nb_kernel_sse2_single_h
#define nb_kernel_sse2_single_h

#include "../nb_kernel.h"


/* List of kernels for this architecture with metadata about them */
extern nb_kernel_info_t
    kernellist_sse2_single[];

/* Length of kernellist_c */
extern int
    kernellist_sse2_single_size;

#endif
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
