/*
<<<<<<< HEAD
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
=======
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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

#ifndef _x2top_h
#define _x2top_h

<<<<<<< HEAD
	
#include <stdio.h>
	
typedef struct {
  char   *elem,*type;
  double q,m;
  int    nbonds;
  char   **bond;
  double *blen;
} t_nm2type;

extern t_nm2type *rd_nm2type(const char *ffdir,int *nnm);
/* Read the name 2 type database. nnm is the number of entries 
 * ff is the force field.
 */

extern void dump_nm2type(FILE *fp,int nnm,t_nm2type nm2t[]);
/* Dump the database for debugging. Can be reread by the program */

extern int nm2type(int nnm,t_nm2type nm2t[],t_symtab *tab,t_atoms *atoms,
		   gpp_atomtype_t atype,int *nbonds,t_params *bond);
/* Try to determine the atomtype (force field dependent) for the atoms 
 * with help of the bond list 
=======

#include <stdio.h>
#include "visibility.h"
typedef struct {
    char    *elem, *type;
    double   q, m;
    int      nbonds;
    char   **bond;
    double  *blen;
} t_nm2type;

GMX_LIBGMXPREPROCESS_EXPORT
extern t_nm2type *rd_nm2type(const char *ffdir, int *nnm);
/* Read the name 2 type database. nnm is the number of entries
 * ff is the force field.
 */

GMX_LIBGMXPREPROCESS_EXPORT
extern void dump_nm2type(FILE *fp, int nnm, t_nm2type nm2t[]);
/* Dump the database for debugging. Can be reread by the program */

GMX_LIBGMXPREPROCESS_EXPORT
extern int nm2type(int nnm, t_nm2type nm2t[], t_symtab *tab, t_atoms *atoms,
                   gpp_atomtype_t atype, int *nbonds, t_params *bond);
/* Try to determine the atomtype (force field dependent) for the atoms
 * with help of the bond list
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */

#endif
