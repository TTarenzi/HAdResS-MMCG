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

#ifndef _gpp_atomtype_h
#define _gpp_atomtype_h

#include <stdio.h>
<<<<<<< HEAD
=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"
#include "macros.h"
#include "grompp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gpp_atomtype *gpp_atomtype_t;

<<<<<<< HEAD
int get_atomtype_type(const char *str,gpp_atomtype_t at);
/* Return atomtype corresponding to case-insensitive str
   or NOTSET if not found */

int get_atomtype_ntypes(gpp_atomtype_t at);
/* Return number of atomtypes */

char *get_atomtype_name(int nt,gpp_atomtype_t at);
/* Return name corresponding to atomtype nt, or NULL if not found */

real get_atomtype_massA(int nt,gpp_atomtype_t at);
real get_atomtype_massB(int nt,gpp_atomtype_t at);
real get_atomtype_qA(int nt,gpp_atomtype_t at);
real get_atomtype_qB(int nt,gpp_atomtype_t at);
real get_atomtype_radius(int nt,gpp_atomtype_t at);
real get_atomtype_vol(int nt,gpp_atomtype_t at);
real get_atomtype_surftens(int nt,gpp_atomtype_t at);
real get_atomtype_gb_radius(int nt,gpp_atomtype_t at);
real get_atomtype_S_hct(int nt,gpp_atomtype_t at);
int get_atomtype_ptype(int nt,gpp_atomtype_t at);
int get_atomtype_batype(int nt,gpp_atomtype_t at);
int get_atomtype_atomnumber(int nt,gpp_atomtype_t at);

/* Return the above variable for atomtype nt, or NOTSET if not found */

real get_atomtype_nbparam(int nt,int param,gpp_atomtype_t at);
/* Similar to the previous but returns the paramth parameter or NOTSET */

=======
GMX_LIBGMXPREPROCESS_EXPORT
int get_atomtype_type(const char *str, gpp_atomtype_t at);
/* Return atomtype corresponding to case-insensitive str
   or NOTSET if not found */

GMX_LIBGMXPREPROCESS_EXPORT
int get_atomtype_ntypes(gpp_atomtype_t at);
/* Return number of atomtypes */

GMX_LIBGMXPREPROCESS_EXPORT
char *get_atomtype_name(int nt, gpp_atomtype_t at);
/* Return name corresponding to atomtype nt, or NULL if not found */

real get_atomtype_massA(int nt, gpp_atomtype_t at);
real get_atomtype_massB(int nt, gpp_atomtype_t at);
real get_atomtype_qA(int nt, gpp_atomtype_t at);
real get_atomtype_qB(int nt, gpp_atomtype_t at);
GMX_LIBGMXPREPROCESS_EXPORT
real get_atomtype_radius(int nt, gpp_atomtype_t at);
GMX_LIBGMXPREPROCESS_EXPORT
real get_atomtype_vol(int nt, gpp_atomtype_t at);
GMX_LIBGMXPREPROCESS_EXPORT
real get_atomtype_surftens(int nt, gpp_atomtype_t at);
GMX_LIBGMXPREPROCESS_EXPORT
real get_atomtype_gb_radius(int nt, gpp_atomtype_t at);
GMX_LIBGMXPREPROCESS_EXPORT
real get_atomtype_S_hct(int nt, gpp_atomtype_t at);
int get_atomtype_ptype(int nt, gpp_atomtype_t at);
int get_atomtype_batype(int nt, gpp_atomtype_t at);
GMX_LIBGMXPREPROCESS_EXPORT
int get_atomtype_atomnumber(int nt, gpp_atomtype_t at);

/* Return the above variable for atomtype nt, or NOTSET if not found */

real get_atomtype_nbparam(int nt, int param, gpp_atomtype_t at);
/* Similar to the previous but returns the paramth parameter or NOTSET */

GMX_LIBGMXPREPROCESS_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
gpp_atomtype_t init_atomtype(void);
/* Return a new atomtype structure */

void done_atomtype(gpp_atomtype_t at);
/* Free the memory in the structure */

<<<<<<< HEAD
int set_atomtype(int nt,gpp_atomtype_t at,t_symtab *tab,
			t_atom *a,const char *name,t_param *nb,
			int bondatomtype,
			real radius,real vol,real surftens,int atomnumber,
			real gb_radius, real S_hct);
/* Set the values of an existing atom type nt. Returns nt on success or
   NOTSET on error. */	

int
set_atomtype_gbparam(gpp_atomtype_t at, int i,
		     real radius,real vol,real surftens,
		     real gb_radius, real S_hct);

int add_atomtype(gpp_atomtype_t at,t_symtab *tab,
			t_atom *a,const char *name,t_param *nb,
			int bondatomtype,
			real radius,real vol,real surftens,real atomnumber,
			real gb_radius, real S_hct);
=======
int set_atomtype(int nt, gpp_atomtype_t at, t_symtab *tab,
                 t_atom *a, const char *name, t_param *nb,
                 int bondatomtype,
                 real radius, real vol, real surftens, int atomnumber,
                 real gb_radius, real S_hct);
/* Set the values of an existing atom type nt. Returns nt on success or
   NOTSET on error. */

int
set_atomtype_gbparam(gpp_atomtype_t at, int i,
                     real radius, real vol, real surftens,
                     real gb_radius, real S_hct);

int add_atomtype(gpp_atomtype_t at, t_symtab *tab,
                 t_atom *a, const char *name, t_param *nb,
                 int bondatomtype,
                 real radius, real vol, real surftens, real atomnumber,
                 real gb_radius, real S_hct);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Add a complete new atom type to an existing atomtype structure. Returns
   the number of the atom type. */

void print_at (FILE * out, gpp_atomtype_t at);
/* Print an atomtype record to a text file */

<<<<<<< HEAD
void renum_atype(t_params plist[],gmx_mtop_t *mtop,
			int *wall_atomtype,
			gpp_atomtype_t at,gmx_bool bVerbose);
			
void copy_atomtype_atomtypes(gpp_atomtype_t atype,t_atomtypes *atypes);
=======
GMX_LIBGMXPREPROCESS_EXPORT
void renum_atype(t_params plist[], gmx_mtop_t *mtop,
                 int *wall_atomtype,
                 gpp_atomtype_t at, gmx_bool bVerbose);

GMX_LIBGMXPREPROCESS_EXPORT
void copy_atomtype_atomtypes(gpp_atomtype_t atype, t_atomtypes *atypes);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Copy from one structure to another */

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _gpp_atomtype_h */








=======
#endif  /* _gpp_atomtype_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2