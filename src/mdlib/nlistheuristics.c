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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
=======
/*
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

<<<<<<< HEAD
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _isnan() */
#include <float.h>
#endif

#include "typedefs.h"
#include "gmx_fatal.h"
#include "vec.h"

void reset_nlistheuristics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    nlh->lt_runav  = 0;
    nlh->lt_runav2 = 0;
=======
#include "typedefs.h"
#include "types/nlistheuristics.h"
#include "gmx_fatal.h"
#include "vec.h"

void reset_nlistheuristics(gmx_nlheur_t *nlh, gmx_large_int_t step)
{
    nlh->lt_runav     = 0;
    nlh->lt_runav2    = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    nlh->step_nscheck = step;
}

void init_nlistheuristics(gmx_nlheur_t *nlh,
<<<<<<< HEAD
                          gmx_bool bGStatEveryStep,gmx_large_int_t step)
{
    nlh->bGStatEveryStep = bGStatEveryStep;
    nlh->nns       = 0;
    nlh->nabnsb    = 0;
    nlh->s1        = 0;
    nlh->s2        = 0;
    nlh->ab        = 0;

    reset_nlistheuristics(nlh,step);
}

void update_nliststatistics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    gmx_large_int_t nl_lt;
    char sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];
=======
                          gmx_bool bGStatEveryStep, gmx_large_int_t step)
{
    nlh->bGStatEveryStep = bGStatEveryStep;
    nlh->nns             = 0;
    nlh->nabnsb          = 0;
    nlh->s1              = 0;
    nlh->s2              = 0;
    nlh->ab              = 0;

    reset_nlistheuristics(nlh, step);
}

void update_nliststatistics(gmx_nlheur_t *nlh, gmx_large_int_t step)
{
    gmx_large_int_t nl_lt;
    char            sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Determine the neighbor list life time */
    nl_lt = step - nlh->step_ns;
    if (debug)
    {
<<<<<<< HEAD
        fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %s steps\n",nlh->nabnsb,gmx_step_str(nl_lt,sbuf));
=======
        fprintf(debug, "%d atoms beyond ns buffer, updating neighbor list after %s steps\n", nlh->nabnsb, gmx_step_str(nl_lt, sbuf));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    nlh->nns++;
    nlh->s1 += nl_lt;
    nlh->s2 += nl_lt*nl_lt;
    nlh->ab += nlh->nabnsb;
    if (nlh->lt_runav == 0)
    {
        nlh->lt_runav  = nl_lt;
        /* Initialize the fluctuation average
         * such that at startup we check after 0 steps.
         */
        nlh->lt_runav2 = sqr(nl_lt/2.0);
    }
    /* Running average with 0.9 gives an exp. history of 9.5 */
    nlh->lt_runav2 = 0.9*nlh->lt_runav2 + 0.1*sqr(nlh->lt_runav - nl_lt);
    nlh->lt_runav  = 0.9*nlh->lt_runav  + 0.1*nl_lt;
    if (nlh->bGStatEveryStep)
    {
        /* Always check the nlist validity */
        nlh->step_nscheck = step;
    }
    else
    {
        /* We check after:  <life time> - 2*sigma
         * The factor 2 is quite conservative,
         * but we assume that with nstlist=-1 the user
         * prefers exact integration over performance.
         */
        nlh->step_nscheck = step
<<<<<<< HEAD
                  + (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)) - 1;
    }
    if (debug)
    {
        fprintf(debug,"nlist life time %s run av. %4.1f sig %3.1f check %s check with -gcom %d\n",
                gmx_step_str(nl_lt,sbuf),nlh->lt_runav,sqrt(nlh->lt_runav2),
                gmx_step_str(nlh->step_nscheck-step+1,sbuf2),
=======
            + (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)) - 1;
    }
    if (debug)
    {
        fprintf(debug, "nlist life time %s run av. %4.1f sig %3.1f check %s check with -gcom %d\n",
                gmx_step_str(nl_lt, sbuf), nlh->lt_runav, sqrt(nlh->lt_runav2),
                gmx_step_str(nlh->step_nscheck-step+1, sbuf2),
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)));
    }
}

<<<<<<< HEAD
void set_nlistheuristics(gmx_nlheur_t *nlh,gmx_bool bReset,gmx_large_int_t step)
=======
void set_nlistheuristics(gmx_nlheur_t *nlh, gmx_bool bReset, gmx_large_int_t step)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    int d;

    if (bReset)
    {
<<<<<<< HEAD
        reset_nlistheuristics(nlh,step);
    }
    else
    {
        update_nliststatistics(nlh,step);
=======
        reset_nlistheuristics(nlh, step);
    }
    else
    {
        update_nliststatistics(nlh, step);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    nlh->step_ns = step;
    /* Initialize the cumulative coordinate scaling matrix */
    clear_mat(nlh->scale_tot);
<<<<<<< HEAD
    for(d=0; d<DIM; d++)
=======
    for (d = 0; d < DIM; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        nlh->scale_tot[d][d] = 1.0;
    }
}
