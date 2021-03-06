/*
<<<<<<< HEAD
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
 * GRoups of Organic Molecules in ACtion for Science
=======
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
#ifndef _iteratedconstraints_h
#define _iteratedconstraints_h

<<<<<<< HEAD
=======
#include "../visibility.h"

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef __cplusplus
extern "C" {
#endif

#if 0
}
/* Hack to make automatic indenting work */
#endif

/* Definitions for convergence of iterated constraints */

/* iterate constraints up to 50 times  */
#define MAXITERCONST       50

/* data type */
typedef struct
{
<<<<<<< HEAD
    real f,fprev,x,xprev;  
    int iter_i;
    gmx_bool bIterate;
    real allrelerr[MAXITERCONST+2];
    int num_close; /* number of "close" violations, caused by limited precision. */
} gmx_iterate_t;

void gmx_iterate_init(gmx_iterate_t *iterate,gmx_bool bIterate);

gmx_bool done_iterating(const t_commrec *cr,FILE *fplog, int nsteps, gmx_iterate_t *iterate, gmx_bool bFirstIterate, real fom, real *newf);
=======
    real     f, fprev, x, xprev;
    int      iter_i;
    gmx_bool bIterationActive;
    real     allrelerr[MAXITERCONST+2];
    int      num_close; /* number of "close" violations, caused by limited precision. */
} gmx_iterate_t;

GMX_LIBMD_EXPORT
void gmx_iterate_init(gmx_iterate_t *iterate, gmx_bool bIterate);

GMX_LIBMD_EXPORT
gmx_bool done_iterating(const t_commrec *cr, FILE *fplog, int nsteps, gmx_iterate_t *iterate, gmx_bool bFirstIterate, real fom, real *newf);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
}
#endif

#endif
