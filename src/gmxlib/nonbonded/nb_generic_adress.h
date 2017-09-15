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
 *                        VERSION 4.0.5
 * Written by Christoph Junghans, Brad Lambeth, and possibly others.
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * All rights reserved.
 
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * Copyright (c) 2011 Christoph Junghans, Sebastian Fritsch
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

#ifndef _nb_generic_adress_h_
#define _nb_generic_adress_h_

<<<<<<< HEAD
=======
#include "nb_kernel.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "types/simple.h"
#include "typedefs.h"

void
<<<<<<< HEAD
gmx_nb_generic_adress_kernel(t_nblist *           nlist,
			 t_forcerec *         fr,
			 t_mdatoms *          mdatoms,
			 real *               x,
			 real *               f,
			 real *               fshift,
			 real *               Vc,
			 real *               Vvdw,
			 real                 tabscale,  
			 real *               VFtab,
			 int *                outeriter,
			 int *                inneriter,
                         gmx_bool                 bCG);

#endif

=======
gmx_nb_generic_adress_kernel(t_nblist *                nlist,
                             rvec *                    xx,
                             rvec *                    ff,
                             t_forcerec *              fr,
                             t_mdatoms *               mdatoms,
                             nb_kernel_data_t *        kernel_data,
                             t_nrnb *                  nrnb);

#endif
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
