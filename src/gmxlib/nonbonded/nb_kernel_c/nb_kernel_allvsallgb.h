<<<<<<< HEAD
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
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifndef _NB_KERNEL_ALLVSALLGB_H
#define _NB_KERNEL_ALLVSALLGB_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "typedefs.h"
<<<<<<< HEAD

void
nb_kernel_allvsallgb(t_forcerec *           fr,
                     t_mdatoms *            mdatoms,
                     t_blocka *             excl,    
                     real *                 x,
                     real *                 f,
                     real *                 Vc,
                     real *                 Vvdw,
                     real *                 vpol,
                     int *                  outeriter,
                     int *                  inneriter,
                     void *                 work);
=======
#include "../nb_kernel.h"

void
nb_kernel_allvsallgb(t_nblist *                nlist,
                     rvec *                    x,
                     rvec *                    f,
                     t_forcerec *              fr,
                     t_mdatoms *               mdatoms,
                     nb_kernel_data_t *        kernel_data,
                     t_nrnb *                  nrnb);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#endif
