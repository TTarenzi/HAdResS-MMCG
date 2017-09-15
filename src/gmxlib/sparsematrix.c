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
 * Green Red Orange Magenta Azure Cyan Skyblue
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "smalloc.h"

#include "sparsematrix.h"




gmx_sparsematrix_t *
gmx_sparsematrix_init(int                    nrow)
{
<<<<<<< HEAD
    int i;
    gmx_sparsematrix_t *A;
    
    snew(A,1);
    
    A->nrow=nrow;
    snew(A->ndata,nrow);
    snew(A->nalloc,nrow);
    snew(A->data,nrow);
        
    for(i=0;i<nrow;i++) 
=======
    int                 i;
    gmx_sparsematrix_t *A;

    snew(A, 1);

    A->nrow = nrow;
    snew(A->ndata, nrow);
    snew(A->nalloc, nrow);
    snew(A->data, nrow);

    for (i = 0; i < nrow; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        A->ndata[i]    = 0;
        A->nalloc[i]   = 0;
        A->data[i]     = NULL;
    }
<<<<<<< HEAD
    return A;    
=======
    return A;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}



void
gmx_sparsematrix_destroy(gmx_sparsematrix_t *   A)
{
    int i;
<<<<<<< HEAD
    
    /* Release each row */
    for(i=0;i<A->nrow;i++) 
    {
        if(A->data[i]!=NULL)
            sfree(A->data[i]);
=======

    /* Release each row */
    for (i = 0; i < A->nrow; i++)
    {
        if (A->data[i] != NULL)
        {
            sfree(A->data[i]);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    /* Release the rowdata arrays */
    sfree(A->ndata);
    sfree(A->nalloc);
    sfree(A->data);
    /* Release matrix structure itself */
<<<<<<< HEAD
    sfree(A);    
=======
    sfree(A);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}



void
gmx_sparsematrix_print(FILE *                 stream,
                       gmx_sparsematrix_t *   A)
{
<<<<<<< HEAD
    int i,j,k;
    
    for(i=0;i<A->nrow;i++)
    {
        if(A->ndata[i]==0) 
        {
            for(j=0;j<A->nrow;j++)
                fprintf(stream," %6.3f",0.0);
        }
        else
        {
            k=0;
            j=0;
            for(j=0;j<A->ndata[i];j++) 
            {
                while(k++<A->data[i][j].col) 
                    fprintf(stream," %6.3f",0.0);
                fprintf(stream," %6.3f",A->data[i][j].value);
            }
            while(k++<A->nrow)
                fprintf(stream," %6.3f",0.0);
        }
        fprintf(stream,"\n");
    }
    
=======
    int i, j, k;

    for (i = 0; i < A->nrow; i++)
    {
        if (A->ndata[i] == 0)
        {
            for (j = 0; j < A->nrow; j++)
            {
                fprintf(stream, " %6.3f", 0.0);
            }
        }
        else
        {
            k = 0;
            j = 0;
            for (j = 0; j < A->ndata[i]; j++)
            {
                while (k++ < A->data[i][j].col)
                {
                    fprintf(stream, " %6.3f", 0.0);
                }
                fprintf(stream, " %6.3f", A->data[i][j].value);
            }
            while (k++ < A->nrow)
            {
                fprintf(stream, " %6.3f", 0.0);
            }
        }
        fprintf(stream, "\n");
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


real
gmx_sparsematrix_value(gmx_sparsematrix_t *    A,
<<<<<<< HEAD
                       int                     row, 
                       int                     col)
{
    gmx_bool found  = FALSE;
    int  i;
    real value;
    
    assert(row<A->nrow);

    value = 0;
    
    /* Find previous value */
    for(i=0;i<A->ndata[row] && (found==FALSE);i++)
    {
        if(A->data[row][i].col==col) 
=======
                       int                     row,
                       int                     col)
{
    gmx_bool found  = FALSE;
    int      i;
    real     value;

    assert(row < A->nrow);

    value = 0;

    /* Find previous value */
    for (i = 0; i < A->ndata[row] && (found == FALSE); i++)
    {
        if (A->data[row][i].col == col)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            found = TRUE;
            value = A->data[row][i].value;
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* value=0 if we didn't find any match */
    return value;
}



void
gmx_sparsematrix_increment_value(gmx_sparsematrix_t *    A,
<<<<<<< HEAD
                                 int                     row, 
=======
                                 int                     row,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                 int                     col,
                                 real                    difference)
{
    gmx_bool found  = FALSE;
<<<<<<< HEAD
    int i;
    
    assert(row<A->nrow);
    
    /* Try to find a previous entry with this row/col */
    for(i=0;i<A->ndata[row] && !found;i++)
    {
        if(A->data[row][i].col==col) 
        {
            found = TRUE;
            A->data[row][i].value += difference;
        }
    }
    
    /* Add a new entry if nothing was found */
    if(!found) 
    {
        /* add the value at the end of the row */
        if(A->ndata[row] == A->nalloc[row]) 
        {
            A->nalloc[row] += 100;
            if(A->data[row]==NULL)
            {
                snew(A->data[row],A->nalloc[row]);
            }
            else
            {
                srenew(A->data[row],A->nalloc[row]);
=======
    int      i;

    assert(row < A->nrow);

    /* Try to find a previous entry with this row/col */
    for (i = 0; i < A->ndata[row] && !found; i++)
    {
        if (A->data[row][i].col == col)
        {
            found                  = TRUE;
            A->data[row][i].value += difference;
        }
    }

    /* Add a new entry if nothing was found */
    if (!found)
    {
        /* add the value at the end of the row */
        if (A->ndata[row] == A->nalloc[row])
        {
            A->nalloc[row] += 100;
            if (A->data[row] == NULL)
            {
                snew(A->data[row], A->nalloc[row]);
            }
            else
            {
                srenew(A->data[row], A->nalloc[row]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
        A->data[row][A->ndata[row]].col   = col;
        /* Previous value was 0.0 */
        A->data[row][A->ndata[row]].value = difference;
        A->ndata[row]++;
<<<<<<< HEAD
    }     
=======
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* Routine to compare column values of two entries, used for quicksort of each row.
 *
 * The data entries to compare are of the type gmx_sparsematrix_entry_t, but quicksort
 * uses void pointers as arguments, so we cast them back internally.
 */
static int
compare_columns(const void *v1, const void *v2)
{
    int c1 = ((gmx_sparsematrix_entry_t *)v1)->col;
    int c2 = ((gmx_sparsematrix_entry_t *)v2)->col;
<<<<<<< HEAD
    
    if(c1<c2)
        return -1;
    else if(c1>c2)
        return 1;
    else 
        return 0;
=======

    if (c1 < c2)
    {
        return -1;
    }
    else if (c1 > c2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


void
gmx_sparsematrix_compress(gmx_sparsematrix_t *    A)
{
<<<<<<< HEAD
    int i,j;
    
    for (i=0;i<A->nrow;i++) 
    {        
        /* Remove last value on this row while it is zero */
        while(A->ndata[i]>0 && A->data[i][A->ndata[i]-1].value==0)
            A->ndata[i]--;
        
        /* Go through values on this row and look for more zero elements */
        for(j=0;j<A->ndata[i];j++)
=======
    int i, j;

    for (i = 0; i < A->nrow; i++)
    {
        /* Remove last value on this row while it is zero */
        while (A->ndata[i] > 0 && A->data[i][A->ndata[i]-1].value == 0)
        {
            A->ndata[i]--;
        }

        /* Go through values on this row and look for more zero elements */
        for (j = 0; j < A->ndata[i]; j++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* If this element was zero, exchange it with the last non-zero
             * element on the row (yes, this will invalidate the sort order)
             */
<<<<<<< HEAD
            if(A->data[i][j].value==0)
=======
            if (A->data[i][j].value == 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                A->data[i][j].value = A->data[i][A->ndata[i]-1].value;
                A->data[i][j].col   = A->data[i][A->ndata[i]-1].col;
                A->ndata[i]--;
            }
        }
        /* Only non-zero elements remaining on this row. Sort them after column index */
        qsort((void *)(A->data[i]),
              A->ndata[i],
              sizeof(gmx_sparsematrix_entry_t),
              compare_columns);
<<<<<<< HEAD
    }        
=======
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


void
gmx_sparsematrix_vector_multiply(gmx_sparsematrix_t *    A,
                                 real *                  x,
                                 real *                  y)
{
<<<<<<< HEAD
    real                        s,v,xi;
    int                         i,j,k;
    gmx_sparsematrix_entry_t *  data; /* pointer to simplify data access */
    
    for (i = 0; i < A->nrow; i ++)
        y[i] = 0;
    
    if(A->compressed_symmetric)
    {
        for (i = 0; i < A->nrow; i ++)
        {
            xi = x[i];
            s = 0.0;
            data = A->data[i];
            
            for (k=0;k<A->ndata[i];k++)
            {
                j = data[k].col;
                v = data[k].value;
                s += v * x[j];
                if(i!=j)
                    y[j] += v * xi; 
            }
            y[i] += s; 
        }    
=======
    real                        s, v, xi;
    int                         i, j, k;
    gmx_sparsematrix_entry_t *  data; /* pointer to simplify data access */

    for (i = 0; i < A->nrow; i++)
    {
        y[i] = 0;
    }

    if (A->compressed_symmetric)
    {
        for (i = 0; i < A->nrow; i++)
        {
            xi   = x[i];
            s    = 0.0;
            data = A->data[i];

            for (k = 0; k < A->ndata[i]; k++)
            {
                j  = data[k].col;
                v  = data[k].value;
                s += v * x[j];
                if (i != j)
                {
                    y[j] += v * xi;
                }
            }
            y[i] += s;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        /* not compressed symmetric storage */
<<<<<<< HEAD
        for (i = 0; i < A->nrow; i ++)
        {
            xi = x[i];
            s = 0.0;
            data = A->data[i];
            
            for (k=0;k<A->ndata[i];k++) 
            {
                j = data[k].col;
                v = data[k].value;
                s += v * x[j];
            }
            y[i] += s; 
        }    
    }
}



=======
        for (i = 0; i < A->nrow; i++)
        {
            xi   = x[i];
            s    = 0.0;
            data = A->data[i];

            for (k = 0; k < A->ndata[i]; k++)
            {
                j  = data[k].col;
                v  = data[k].value;
                s += v * x[j];
            }
            y[i] += s;
        }
    }
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
