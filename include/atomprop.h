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
 * Gromacs Runs On Most of All Computer Systems
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

#ifndef _atomprop_h
#define _atomprop_h
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef __cplusplus
extern "C" {
#endif

#include "index.h"

/* Abstract type for the atom property database */
typedef struct gmx_atomprop *gmx_atomprop_t;

<<<<<<< HEAD
enum { epropMass, epropVDW, epropDGsol, epropElectroneg, epropElement, 
       epropNR };

gmx_atomprop_t gmx_atomprop_init(void);
/* Initializes and returns the atom properties struct */

void gmx_atomprop_destroy(gmx_atomprop_t aps);
/* Get rid of memory after use */

char *gmx_atomprop_element(gmx_atomprop_t aps,int atomnumber);

int gmx_atomprop_atomnumber(gmx_atomprop_t aps,const char *element);

gmx_bool gmx_atomprop_query(gmx_atomprop_t aps,
                        int eprop,const char *resnm,const char *atomnm,
                        real *value);
=======
enum {
    epropMass, epropVDW, epropDGsol, epropElectroneg, epropElement,
    epropNR
};

GMX_LIBGMX_EXPORT
gmx_atomprop_t gmx_atomprop_init(void);
/* Initializes and returns the atom properties struct */

GMX_LIBGMX_EXPORT
void gmx_atomprop_destroy(gmx_atomprop_t aps);
/* Get rid of memory after use */

char *gmx_atomprop_element(gmx_atomprop_t aps, int atomnumber);

int gmx_atomprop_atomnumber(gmx_atomprop_t aps, const char *element);

GMX_LIBGMX_EXPORT
gmx_bool gmx_atomprop_query(gmx_atomprop_t aps,
                            int eprop, const char *resnm, const char *atomnm,
                            real *value);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Extract a value from the database. Returns TRUE on succes,
 * FALSE otherwise. In the latter case, value is a deafult value.
 * The first time this function is called for this property
 * the database will be read.
 */

#ifdef __cplusplus
}
#endif


#endif