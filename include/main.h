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

#ifndef _main_h
#define _main_h


#include <stdio.h>
<<<<<<< HEAD
=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "network.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
gmx_bool gmx_parallel_env_initialized(void); 
/* 1 when running in a parallel environment, so could also be 1 if
   mdrun was started with: mpirun -np 1.
  
   Use this function only to check whether a parallel environment has
   been initialized, for example when checking whether gmx_finalize()
   needs to be called. Use PAR(cr) to check whether the simulation actually
   has more than one node/thread. */


void gmx_log_open(const char *fn,const t_commrec *cr,
                          gmx_bool bMasterOnly, unsigned long Flags, FILE**);
=======
char *gmx_gethostname(char *name, size_t len);
/* Sets the hostname to the value given by gethostname, if available,
 * and to "unknown" otherwise. name should have at least size len.
 * Returns name.
 */

GMX_LIBGMX_EXPORT
void gmx_log_open(const char *fn, const t_commrec *cr,
                  gmx_bool bMasterOnly, gmx_bool bAppendFiles, FILE**);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Open the log file, if necessary (nprocs > 1) the logfile name is
 * communicated around the ring.
 */

<<<<<<< HEAD
void gmx_log_close(FILE *fp);
/* Close the log file */

void check_multi_int(FILE *log,const gmx_multisim_t *ms,
			    int val,const char *name);
void check_multi_large_int(FILE *log,const gmx_multisim_t *ms,
                           gmx_large_int_t val,const char *name);
/* Check if val is the same on all processors for a mdrun -multi run
 * The string name is used to print to the log file and in a fatal error
 * if the val's don't match.
 */

=======
GMX_LIBGMX_EXPORT
void gmx_log_close(FILE *fp);
/* Close the log file */

GMX_LIBGMX_EXPORT
void check_multi_int(FILE *log, const gmx_multisim_t *ms,
                     int val, const char *name,
                     gmx_bool bQuiet);
GMX_LIBGMX_EXPORT
void check_multi_large_int(FILE *log, const gmx_multisim_t *ms,
                           gmx_large_int_t val, const char *name,
                           gmx_bool bQuiet);
/* Check if val is the same on all processors for a mdrun -multi run
 * The string name is used to print to the log file and in a fatal error
 * if the val's don't match. If bQuiet is true and the check passes,
 * no output is written.
 */

GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
void init_multisystem(t_commrec *cr, int nsim, char **multidirs,
                      int nfile, const t_filenm fnm[], gmx_bool bParFn);
/* Splits the communication into nsim separate simulations
 * and creates a communication structure between the master
 * these simulations.
 * If bParFn is set, the nodeid is appended to the tpx and each output file.
 */

<<<<<<< HEAD
t_commrec *init_par(int *argc,char ***argv_ptr);
=======
GMX_LIBGMX_EXPORT
t_commrec *init_par(int *argc, char ***argv_ptr);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Initiate the parallel computer. Return the communication record
 * (see network.h). The command line arguments are communicated so that they can be
 * parsed on each processor.
 * Arguments are the number of command line arguments, and a pointer to the
<<<<<<< HEAD
 * array of argument strings.
 */

t_commrec *init_par_threads(const t_commrec *cro);
/* Initialize communication records for thread-parallel simulations. 
   Must be called on all threads before any communication takes place by 
   the individual threads. Copies the original commrec to 
   thread-local versions (a small memory leak results because we don't 
   deallocate the old shared version).  */

t_commrec *init_cr_nopar(void);
/* Returns t_commrec for non-parallel functionality */

=======
 * array of argument strings. Both are allowed to be NULL.
 */

GMX_LIBGMX_EXPORT
t_commrec *init_par_threads(const t_commrec *cro);
/* Initialize communication records for thread-parallel simulations.
   Must be called on all threads before any communication takes place by
   the individual threads. Copies the original commrec to
   thread-local versions (a small memory leak results because we don't
   deallocate the old shared version).  */

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _main_h */
=======
#endif  /* _main_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
