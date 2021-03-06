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

#ifndef _confio_h
#define _confio_h

<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"

/* For reading coordinate files it is assumed that enough memory
 * has been allocated beforehand.
 */
#ifdef __cplusplus
extern "C" {
#endif
<<<<<<< HEAD
  
void init_t_atoms(t_atoms *atoms, int natoms, gmx_bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */  

int read_g96_conf(FILE *fp,const char *infile,t_trxframe *fr);
/* read a Gromos96 coordinate or trajectory file,                       *
 * returns the number of atoms                                          *
 * sets what's in the frame in info                                     *  
 * read from fp, infile is only needed for error messages               *   
 * nwanted is the number of wanted coordinates,                         *
 * set this to -1 if you want to know the number of atoms in the file   *
 * title, atoms, x, v can all be NULL, in which case they won't be read */

void write_g96_conf(FILE *out,t_trxframe *fr,int nindex,atom_id *index);
/* write a Gromos96 coordinate file or trajectory frame *
 * index can be NULL                                    */

gmx_bool gro_next_x_or_v(FILE *status,t_trxframe *fr);
int gro_first_x_or_v(FILE *status,t_trxframe *fr);
/* read first/next x and/or v frame from gro file */

void write_hconf_indexed_p(FILE *out,const char *title,t_atoms *atoms,
				  int nx,atom_id index[],int ndec,
				  rvec *x,rvec *v,matrix box);
		
void write_hconf_p(FILE *out,const char *title,t_atoms *atoms, int ndec,
			  rvec *x,rvec *v,matrix box); 
/* Write a Gromos file with precision ndec: number of decimal places in x,
 * v has one place more. */ 

void write_sto_conf_indexed(const char *outfile,const char *title,
			    t_atoms *atoms, 
			    rvec x[],rvec *v,int ePBC,matrix box,
			    atom_id nindex,atom_id index[]);
/* like write_sto_conf, but indexed */ 

void write_sto_conf(const char *outfile,const char *title,
			   t_atoms *atoms, 
			   rvec x[],rvec *v,int ePBC,matrix box);
/* write atoms, x, v (if .gro and not NULL) and box (if not NULL) 
 * to an STO (.gro or .pdb) file */ 

void write_sto_conf_mtop(const char *outfile,const char *title,
				gmx_mtop_t *mtop,
				rvec x[],rvec *v,int ePBC,matrix box);
/* As write_sto_conf, but uses a gmx_mtop_t struct */

void get_stx_coordnum (const char *infile,int *natoms);
/* read the number of atoms from an STX file */

void read_stx_conf(const char *infile,char *title,
			  t_atoms *atoms, 
			  rvec x[],rvec *v,int *ePBC,matrix box);
=======

int read_g96_conf(FILE *fp, const char *infile, t_trxframe *fr, char *line);
/* read a Gromos96 coordinate or trajectory file,                       *
 * returns the number of atoms                                          *
 * sets what's in the frame in info                                     *
 * read from fp, infile is only needed for error messages               *
 * nwanted is the number of wanted coordinates,                         *
 * set this to -1 if you want to know the number of atoms in the file   *
 * title, atoms, x, v can all be NULL, in which case they won't be read *
 * line holds the previous line for trajectory reading                  */

GMX_LIBGMX_EXPORT
void write_g96_conf(FILE *out, t_trxframe *fr, int nindex, atom_id *index);
/* write a Gromos96 coordinate file or trajectory frame *
 * index can be NULL                                    */

gmx_bool gro_next_x_or_v(FILE *status, t_trxframe *fr);
int gro_first_x_or_v(FILE *status, t_trxframe *fr);
/* read first/next x and/or v frame from gro file */

void write_hconf_indexed_p(FILE *out, const char *title, t_atoms *atoms,
                           int nx, atom_id index[], int ndec,
                           rvec *x, rvec *v, matrix box);

GMX_LIBGMX_EXPORT
void write_hconf_p(FILE *out, const char *title, t_atoms *atoms, int ndec,
                   rvec *x, rvec *v, matrix box);
/* Write a Gromos file with precision ndec: number of decimal places in x,
 * v has one place more. */

GMX_LIBGMX_EXPORT
void write_sto_conf_indexed(const char *outfile, const char *title,
                            t_atoms *atoms,
                            rvec x[], rvec *v, int ePBC, matrix box,
                            atom_id nindex, atom_id index[]);
/* like write_sto_conf, but indexed */

GMX_LIBGMX_EXPORT
void write_sto_conf(const char *outfile, const char *title,
                    t_atoms *atoms,
                    rvec x[], rvec *v, int ePBC, matrix box);
/* write atoms, x, v (if .gro and not NULL) and box (if not NULL)
 * to an STO (.gro or .pdb) file */

GMX_LIBGMX_EXPORT
void write_sto_conf_mtop(const char *outfile, const char *title,
                         gmx_mtop_t *mtop,
                         rvec x[], rvec *v, int ePBC, matrix box);
/* As write_sto_conf, but uses a gmx_mtop_t struct */

GMX_LIBGMX_EXPORT
void get_stx_coordnum (const char *infile, int *natoms);
/* read the number of atoms from an STX file */

GMX_LIBGMX_EXPORT
void read_stx_conf(const char *infile, char *title,
                   t_atoms *atoms,
                   rvec x[], rvec *v, int *ePBC, matrix box);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Read atoms, x, v and box from an STX file.
 * If ePBC!=NULL return the type of pbc in *ePBC or -1 if unknown.
 */

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _confio_h */
=======
#endif  /* _confio_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
