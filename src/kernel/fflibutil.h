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

#ifndef _fflibutil_h
#define _fflibutil_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
<<<<<<< HEAD
=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const char *fflib_forcefield_dir_ext();
/* Returns the name of the force field directory extension */

extern const char *fflib_forcefield_itp();
/* Returns the name of the main forcefield itp file */

extern const char *fflib_forcefield_doc();
/* Returns the name of the forcefield documentation file */

<<<<<<< HEAD
extern void fflib_filename_base(const char *filename,char *filebase,int maxlen);
=======
extern void fflib_filename_base(const char *filename, char *filebase, int maxlen);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Return the base file name of filename in base,
 * i.e. remove path and extension, if present.
 * base should be at least of size maxlen.
 */

<<<<<<< HEAD
extern int fflib_search_file_end(const char *ffdir,
				 const char *file_end,
				 gmx_bool bFatalError,
				 char ***filenames);
=======
GMX_LIBGMXPREPROCESS_EXPORT
extern int fflib_search_file_end(const char *ffdir,
                                 const char *file_end,
                                 gmx_bool    bFatalError,
                                 char     ***filenames);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Search for files ending on file_end in the force field directory fflib.
 * fflib should be in the GROMACS lib.path.
 * Return the number of files and the file names in filenames.
 */

<<<<<<< HEAD
extern int fflib_search_file_in_dirend(const char *filename,const char *dirend,
				       char ***dirnames);
=======
extern int fflib_search_file_in_dirend(const char *filename, const char *dirend,
                                       char ***dirnames);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Search for files with name filename in subdirectories with names
 * ending on dirend.
 * Return the number of files and the directory names in dirnames.
 */
<<<<<<< HEAD
extern gmx_bool fflib_fexist(const char *file);
/* Check if a file exists in the force field library */

=======
GMX_LIBGMXPREPROCESS_EXPORT
extern gmx_bool fflib_fexist(const char *file);
/* Check if a file exists in the force field library */

GMX_LIBGMXPREPROCESS_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
extern FILE *fflib_open(const char *file);
/* Open force field library file "file" for reading.
 * "file" should contain the whole path to the force field library,
 * either absolute or relative to the current dir.
 */

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _fflibutil_h */
=======
#endif  /* _fflibutil_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
