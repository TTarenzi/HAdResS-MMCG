<<<<<<< HEAD
/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
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
 * GROningen Mixture of Alchemy and Childrens' Stories
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

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "mshift.h"
#include "pbc.h"
#include "gstat.h"
#include "futil.h"
<<<<<<< HEAD
#include "vec.h"	

typedef struct {
    int     natoms;
=======
#include "vec.h"

typedef struct {
    int      natoms;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_graph *gr;
} rmpbc_graph_t;

typedef struct gmx_rmpbc {
    t_idef        *idef;
<<<<<<< HEAD
    int           natoms_init;
    int           ePBC;
    int           ngraph;
    rmpbc_graph_t *graph;
} koeiepoep;

static t_graph *gmx_rmpbc_get_graph(gmx_rmpbc_t gpbc,int ePBC,int natoms)
{
    int           i;
=======
    int            natoms_init;
    int            ePBC;
    int            ngraph;
    rmpbc_graph_t *graph;
} koeiepoep;

static t_graph *gmx_rmpbc_get_graph(gmx_rmpbc_t gpbc, int ePBC, int natoms)
{
    int            i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rmpbc_graph_t *gr;

    if (ePBC == epbcNONE
        || NULL == gpbc
        || NULL == gpbc->idef
        || gpbc->idef->ntypes <= 0)
    {
        return NULL;
    }

    gr = NULL;
<<<<<<< HEAD
    for(i=0; i<gpbc->ngraph; i++)
=======
    for (i = 0; i < gpbc->ngraph; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (natoms == gpbc->graph[i].natoms)
        {
            gr = &gpbc->graph[i];
        }
    }
    if (gr == NULL)
    {
        /* We'd like to check with the number of atoms in the topology,
         * but we don't have that available.
         * So we check against the number of atoms that gmx_rmpbc_init
         * was called with.
         */
        if (natoms > gpbc->natoms_init)
        {
<<<<<<< HEAD
            gmx_fatal(FARGS,"Structure or trajectory file has more atoms (%d) than the topology (%d)",natoms,gpbc->natoms_init);
        }
        gpbc->ngraph++;
        srenew(gpbc->graph,gpbc->ngraph);
        gr = &gpbc->graph[gpbc->ngraph-1];
        gr->natoms = natoms;
        gr->gr     = mk_graph(NULL,gpbc->idef,0,natoms,FALSE,FALSE);
=======
            gmx_fatal(FARGS, "Structure or trajectory file has more atoms (%d) than the topology (%d)", natoms, gpbc->natoms_init);
        }
        gpbc->ngraph++;
        srenew(gpbc->graph, gpbc->ngraph);
        gr         = &gpbc->graph[gpbc->ngraph-1];
        gr->natoms = natoms;
        gr->gr     = mk_graph(NULL, gpbc->idef, 0, natoms, FALSE, FALSE);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return gr->gr;
}

<<<<<<< HEAD
gmx_rmpbc_t gmx_rmpbc_init(t_idef *idef,int ePBC,int natoms,
                           matrix box)
{
    gmx_rmpbc_t gpbc;
  
    snew(gpbc,1);

    gpbc->natoms_init = natoms;
  
=======
gmx_rmpbc_t gmx_rmpbc_init(t_idef *idef, int ePBC, int natoms,
                           matrix box)
{
    gmx_rmpbc_t gpbc;

    snew(gpbc, 1);

    gpbc->natoms_init = natoms;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* This sets pbc when we now it,
     * otherwise we guess it from the instantaneous box in the trajectory.
     */
    gpbc->ePBC = ePBC;

    gpbc->idef = idef;
    if (gpbc->idef->ntypes <= 0)
    {
        fprintf(stderr,
                "\n"
<<<<<<< HEAD
                "WARNING: if there are broken molecules in the trajectory file,\n"
                "         they can not be made whole without a run input file\n\n");
=======
                "WARNING: If there are molecules in the input trajectory file\n"
                "         that are broken across periodic boundaries, they\n"
                "         cannot be made whole (or treated as whole) without\n"
                "         you providing a run input file.\n\n");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return gpbc;
}

void gmx_rmpbc_done(gmx_rmpbc_t gpbc)
{
    int i;

    if (NULL != gpbc)
    {
<<<<<<< HEAD
        for(i=0; i<gpbc->ngraph; i++)
=======
        for (i = 0; i < gpbc->ngraph; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            done_graph(gpbc->graph[i].gr);
        }
        if (gpbc->graph != NULL)
        {
            sfree(gpbc->graph);
        }
    }
}

<<<<<<< HEAD
static int gmx_rmpbc_ePBC(gmx_rmpbc_t gpbc,matrix box)
=======
static int gmx_rmpbc_ePBC(gmx_rmpbc_t gpbc, matrix box)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    if (NULL != gpbc && gpbc->ePBC >= 0)
    {
        return gpbc->ePBC;
    }
    else
    {
        return guess_ePBC(box);
    }
}

<<<<<<< HEAD
void gmx_rmpbc(gmx_rmpbc_t gpbc,int natoms,matrix box,rvec x[])
{
    int     ePBC;
    t_graph *gr;
    
    ePBC = gmx_rmpbc_ePBC(gpbc,box);
    gr = gmx_rmpbc_get_graph(gpbc,ePBC,natoms);
    if (gr != NULL)
    {
        mk_mshift(stdout,gr,ePBC,box,x);
        shift_x(gr,box,x,x);
    }
}

void gmx_rmpbc_copy(gmx_rmpbc_t gpbc,int natoms,matrix box,rvec x[],rvec x_s[])
{
    int     ePBC;
    t_graph *gr;
    int     i;

    ePBC = gmx_rmpbc_ePBC(gpbc,box);
    gr = gmx_rmpbc_get_graph(gpbc,ePBC,natoms);
    if (gr != NULL)
    {
        mk_mshift(stdout,gr,ePBC,box,x);
        shift_x(gr,box,x,x_s);
    }
    else
    {
        for(i=0; i<natoms; i++)
        {
            copy_rvec(x[i],x_s[i]);
=======
void gmx_rmpbc(gmx_rmpbc_t gpbc, int natoms, matrix box, rvec x[])
{
    int      ePBC;
    t_graph *gr;

    ePBC = gmx_rmpbc_ePBC(gpbc, box);
    gr   = gmx_rmpbc_get_graph(gpbc, ePBC, natoms);
    if (gr != NULL)
    {
        mk_mshift(stdout, gr, ePBC, box, x);
        shift_self(gr, box, x);
    }
}

void gmx_rmpbc_copy(gmx_rmpbc_t gpbc, int natoms, matrix box, rvec x[], rvec x_s[])
{
    int      ePBC;
    t_graph *gr;
    int      i;

    ePBC = gmx_rmpbc_ePBC(gpbc, box);
    gr   = gmx_rmpbc_get_graph(gpbc, ePBC, natoms);
    if (gr != NULL)
    {
        mk_mshift(stdout, gr, ePBC, box, x);
        shift_x(gr, box, x, x_s);
    }
    else
    {
        for (i = 0; i < natoms; i++)
        {
            copy_rvec(x[i], x_s[i]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
}

<<<<<<< HEAD
void gmx_rmpbc_trxfr(gmx_rmpbc_t gpbc,t_trxframe *fr)
{
    int     ePBC;
=======
void gmx_rmpbc_trxfr(gmx_rmpbc_t gpbc, t_trxframe *fr)
{
    int      ePBC;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_graph *gr;

    if (fr->bX && fr->bBox)
    {
<<<<<<< HEAD
        ePBC = gmx_rmpbc_ePBC(gpbc,fr->box);
        gr = gmx_rmpbc_get_graph(gpbc,ePBC,fr->natoms);
        if (gr != NULL)
        {
            mk_mshift(stdout,gr,ePBC,fr->box,fr->x);
            shift_x(gr,fr->box,fr->x,fr->x);
=======
        ePBC = gmx_rmpbc_ePBC(gpbc, fr->box);
        gr   = gmx_rmpbc_get_graph(gpbc, ePBC, fr->natoms);
        if (gr != NULL)
        {
            mk_mshift(stdout, gr, ePBC, fr->box, fr->x);
            shift_self(gr, fr->box, fr->x);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
}

<<<<<<< HEAD
void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
    real dist;
    int  n,m,d;
  
    /* check periodic boundary */
    for(n=1;(n<atoms->nr);n++)
    {
        for(m=DIM-1; m>=0; m--)
=======
void rm_gropbc(t_atoms *atoms, rvec x[], matrix box)
{
    real dist;
    int  n, m, d;

    /* check periodic boundary */
    for (n = 1; (n < atoms->nr); n++)
    {
        for (m = DIM-1; m >= 0; m--)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            dist = x[n][m]-x[n-1][m];
            if (fabs(dist) > 0.9*box[m][m])
            {
<<<<<<< HEAD
                if ( dist >  0 )
                {
                    for(d=0; d<=m; d++)
=======
                if (dist >  0)
                {
                    for (d = 0; d <= m; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        x[n][d] -= box[m][d];
                    }
                }
                else
                {
<<<<<<< HEAD
                    for(d=0; d<=m; d++)
=======
                    for (d = 0; d <= m; d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        x[n][d] += box[m][d];
                    }
                }
<<<<<<< HEAD
            } 	
        }
    }
}

=======
            }
        }
    }
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2