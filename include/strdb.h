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

#ifndef _strdb_h
#define _strdb_h

<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
gmx_bool get_a_line(FILE *fp,char line[],int n);
/* Read a line of at most n characters form *fp to line. 
 * Comment ';...' and leading spaces are removed, empty lines are skipped.
 * Return FALSE when eof. 
 */

gmx_bool get_header(char line[],char header[]);
=======
GMX_LIBGMX_EXPORT
gmx_bool get_a_line(FILE *fp, char line[], int n);
/* Read a line of at most n characters form *fp to line.
 * Comment ';...' and leading spaces are removed, empty lines are skipped.
 * Return FALSE when eof.
 */

GMX_LIBGMX_EXPORT
gmx_bool get_header(char line[], char header[]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Read a header between '[' and ']' from line to header.
 * Returns FALSE no header is found.
 */

<<<<<<< HEAD
int fget_lines(FILE *in,char ***strings);
=======
int fget_lines(FILE *in, char ***strings);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Read an array of lines from file in. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * return the number of strings.
 */
<<<<<<< HEAD
int get_lines(const char *db,char ***strings);
/* Open file db, or if non-existant file $GMXLIB/db and read strings 
 * return the number of strings.
 */

int search_str(int nstr,char **str,char *key);
=======
GMX_LIBGMX_EXPORT
int get_lines(const char *db, char ***strings);
/* Open file db, or if non-existant file $GMXLIB/db and read strings
 * return the number of strings.
 */

GMX_LIBGMX_EXPORT
int search_str(int nstr, char **str, char *key);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Search an array of strings for key, return the index if found
 * -1 if not found.
 */

<<<<<<< HEAD
int get_strings(const char *db,char ***strings);
=======
GMX_LIBGMX_EXPORT
int get_strings(const char *db, char ***strings);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Read an array of strings from file db or $GMXLIB/db. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * return the number of strings.
 */
<<<<<<< HEAD
int get_file(const char *db,char ***strings);
/* Read an array of strings from file db or $GMXLIB/db. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * Does not need number of lines as first line in the file. 
=======
GMX_LIBGMX_EXPORT
int get_file(const char *db, char ***strings);
/* Read an array of strings from file db or $GMXLIB/db. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * Does not need number of lines as first line in the file.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * return the number of strings.
 */

#ifdef __cplusplus
}
#endif


<<<<<<< HEAD
#endif	/* _strdb_h */
=======
#endif  /* _strdb_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
