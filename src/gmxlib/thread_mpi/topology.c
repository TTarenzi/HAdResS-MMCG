/*
<<<<<<< HEAD
This source code file is part of thread_mpi.  
Written by Sander Pronk, Erik Lindahl, and possibly others. 

Copyright (c) 2009, Sander Pronk, Erik Lindahl.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

If you want to redistribute modifications, please consider that
scientific software is very special. Version control is crucial -
bugs must be traceable. We will be happy to consider code for
inclusion in the official distribution, but derived work should not
be called official thread_mpi. Details are found in the README & COPYING
files.
*/
=======
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2009, Sander Pronk, Erik Lindahl.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   If you want to redistribute modifications, please consider that
   scientific software is very special. Version control is crucial -
   bugs must be traceable. We will be happy to consider code for
   inclusion in the official distribution, but derived work should not
   be called official thread_mpi. Details are found in the README & COPYING
   files.
 */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>


#include "impl.h"


/* topology functions */
int tMPI_Topo_test(tMPI_Comm comm, int *status)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Topo_test(%p, %p)", comm, status);
#endif

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (comm->cart)
<<<<<<< HEAD
        *status=TMPI_CART;
    /*else if (comm->graph)
        status=MPI_GRAPH;*/
    else 
        *status=TMPI_UNDEFINED;
=======
    {
        *status = TMPI_CART;
    }
    /*else if (comm->graph)
        status=MPI_GRAPH;*/
    else
    {
        *status = TMPI_UNDEFINED;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return TMPI_SUCCESS;
}

int tMPI_Cartdim_get(tMPI_Comm comm, int *ndims)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cartdim_get(%p, %p)", comm, ndims);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
<<<<<<< HEAD
    if (!comm->cart || comm->cart->ndims==0)
    {
        return TMPI_SUCCESS;
    }
    *ndims=comm->cart->ndims;
=======
    if (!comm->cart || comm->cart->ndims == 0)
    {
        return TMPI_SUCCESS;
    }
    *ndims = comm->cart->ndims;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return TMPI_SUCCESS;
}


int tMPI_Cart_get(tMPI_Comm comm, int maxdims, int *dims, int *periods,
                  int *coords)
{
    int i;
<<<<<<< HEAD
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_get(%p, %d, %p, %p, %p)", comm, maxdims, 
                       dims, periods, coords);
=======
    int myrank = tMPI_Comm_seek_rank(comm, tMPI_Get_current());

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_get(%p, %d, %p, %p, %p)", comm, maxdims,
                     dims, periods, coords);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
<<<<<<< HEAD
    if (!comm->cart || comm->cart->ndims==0)
        return TMPI_SUCCESS;

    tMPI_Cart_coords(comm, myrank, maxdims, coords);

    for(i=0;i<comm->cart->ndims;i++)
    {
        if (i>=maxdims)
        {
            return tMPI_Error(comm, TMPI_ERR_DIMS);
        }
        dims[i]=comm->cart->dims[i];
        periods[i]=comm->cart->periods[i];
=======
    if (!comm->cart || comm->cart->ndims == 0)
    {
        return TMPI_SUCCESS;
    }

    tMPI_Cart_coords(comm, myrank, maxdims, coords);

    for (i = 0; i < comm->cart->ndims; i++)
    {
        if (i >= maxdims)
        {
            return tMPI_Error(comm, TMPI_ERR_DIMS);
        }
        dims[i]    = comm->cart->dims[i];
        periods[i] = comm->cart->periods[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return TMPI_SUCCESS;
}

int tMPI_Cart_rank(tMPI_Comm comm, int *coords, int *rank)
{
<<<<<<< HEAD
    int i,mul=1,ret=0;
=======
    int i, mul = 1, ret = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_get(%p, %p, %p)", comm, coords, rank);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
<<<<<<< HEAD
    if (!comm->cart || comm->cart->ndims==0)
        return TMPI_SUCCESS;

    /* because of row-major ordering, we count the dimensions down */
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        int rcoord=coords[i];
=======
    if (!comm->cart || comm->cart->ndims == 0)
    {
        return TMPI_SUCCESS;
    }

    /* because of row-major ordering, we count the dimensions down */
    for (i = comm->cart->ndims-1; i >= 0; i--)
    {
        int rcoord = coords[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (comm->cart->periods[i])
        {
            /* apply periodic boundary conditions */
            rcoord = rcoord % comm->cart->dims[i];
            if (rcoord < 0)
<<<<<<< HEAD
                rcoord += comm->cart->dims[i];
=======
            {
                rcoord += comm->cart->dims[i];
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        else
        {
            if (rcoord < 0 || rcoord >= comm->cart->dims[i])
            {
                return tMPI_Error(comm, TMPI_ERR_DIMS);
            }
        }
        ret += mul*rcoord;
        mul *= comm->cart->dims[i];
    }
<<<<<<< HEAD
    *rank=ret;
=======
    *rank = ret;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return TMPI_SUCCESS;
}

int tMPI_Cart_coords(tMPI_Comm comm, int rank, int maxdims, int *coords)
{
    int i;
<<<<<<< HEAD
    int rank_left=rank;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_coords(%p, %d, %d, %p)", comm, rank, maxdims, 
=======
    int rank_left = rank;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_coords(%p, %d, %d, %p)", comm, rank, maxdims,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     coords);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
<<<<<<< HEAD
    if (!comm->cart || comm->cart->ndims==0)
        return TMPI_SUCCESS;
=======
    if (!comm->cart || comm->cart->ndims == 0)
    {
        return TMPI_SUCCESS;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (maxdims < comm->cart->ndims)
    {
        return tMPI_Error(comm, TMPI_ERR_DIMS);
    }

    /* again, row-major ordering */
<<<<<<< HEAD
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        coords[i]=rank_left%comm->cart->dims[i];
        rank_left /= comm->cart->dims[i];
    }   
=======
    for (i = comm->cart->ndims-1; i >= 0; i--)
    {
        coords[i]  = rank_left%comm->cart->dims[i];
        rank_left /= comm->cart->dims[i];
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return TMPI_SUCCESS;
}



<<<<<<< HEAD
int tMPI_Cart_map(tMPI_Comm comm, int ndims, int *dims, int *periods, 
                  int *newrank)
{
    /* this function doesn't actually do anything beyond returning the current 
       rank (or TMPI_UNDEFINED if it doesn't fit in the new topology */
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int Ntot=1;
    int i;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_map(%p, %d, %p, %p, %p)", comm, ndims, dims, 
=======
int tMPI_Cart_map(tMPI_Comm comm, int ndims, int *dims, int *periods,
                  int *newrank)
{
    /* this function doesn't actually do anything beyond returning the current
       rank (or TMPI_UNDEFINED if it doesn't fit in the new topology */
    int myrank = tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int Ntot   = 1;
    int i;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_map(%p, %d, %p, %p, %p)", comm, ndims, dims,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     periods, newrank);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!periods)
    {
        return tMPI_Error(comm, TMPI_ERR_DIMS);
    }
<<<<<<< HEAD
 
    /* calculate the total number of procs in cartesian comm */
    for(i=0;i<ndims;i++)
=======

    /* calculate the total number of procs in cartesian comm */
    for (i = 0; i < ndims; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        Ntot *= dims[i];
    }

    if (myrank >= Ntot)
    {
<<<<<<< HEAD
        *newrank=TMPI_UNDEFINED;
    }
    else
    {
        *newrank=myrank;
=======
        *newrank = TMPI_UNDEFINED;
    }
    else
    {
        *newrank = myrank;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return TMPI_SUCCESS;
}


/* initialize Cartesian topology info in comm. If ndims==0, dims and periods
   are not referenced */
<<<<<<< HEAD
static void tMPI_Cart_init(tMPI_Comm *comm_cart, int ndims, int *dims, 
                           int *periods)
{
    int newrank=-1;
=======
static void tMPI_Cart_init(tMPI_Comm *comm_cart, int ndims, int *dims,
                           int *periods)
{
    int newrank = -1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    int i;

    if (*comm_cart)
    {
        tMPI_Comm_rank(*comm_cart, &newrank);
    }

<<<<<<< HEAD
    if (newrank==0)
    {
        (*comm_cart)->cart=(struct cart_topol*)tMPI_Malloc(
                                            sizeof(struct cart_topol));
        (*comm_cart)->cart->dims=(int*)tMPI_Malloc(ndims*sizeof(int));
        (*comm_cart)->cart->periods=(int*)tMPI_Malloc(ndims*sizeof(int));
        (*comm_cart)->cart->ndims=ndims;
        for(i=0;i<ndims;i++)
        {
            (*comm_cart)->cart->dims[i]=dims[i];
            (*comm_cart)->cart->periods[i]=periods[i];
        }
    }

    /* and we add a barrier to make sure the cart object is seen by 
=======
    if (newrank == 0)
    {
        (*comm_cart)->cart = (struct cart_topol*)tMPI_Malloc(
                    sizeof(struct cart_topol));
        (*comm_cart)->cart->dims    = (int*)tMPI_Malloc(ndims*sizeof(int));
        (*comm_cart)->cart->periods = (int*)tMPI_Malloc(ndims*sizeof(int));
        (*comm_cart)->cart->ndims   = ndims;
        for (i = 0; i < ndims; i++)
        {
            (*comm_cart)->cart->dims[i]    = dims[i];
            (*comm_cart)->cart->periods[i] = periods[i];
        }
    }

    /* and we add a barrier to make sure the cart object is seen by
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
       every thread that is part of the new communicator */
    if (*comm_cart)
    {
        tMPI_Barrier_wait( &( (*comm_cart)->barrier) );
    }
}

void tMPI_Cart_destroy(struct cart_topol *cart)
{
    if (cart)
    {
        free(cart->dims);
        free(cart->periods);
    }
}

int tMPI_Cart_create(tMPI_Comm comm_old, int ndims, int *dims, int *periods,
                     int reorder, tMPI_Comm *comm_cart)
{
<<<<<<< HEAD
    int myrank=tMPI_Comm_seek_rank(comm_old, tMPI_Get_current());
    int key=myrank;
    int color=0;
    int Ntot=1;
    int i;
    

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_create(%p, %d, %p, %p, %d, %p)", comm_old, 
=======
    int myrank = tMPI_Comm_seek_rank(comm_old, tMPI_Get_current());
    int key    = myrank;
    int color  = 0;
    int Ntot   = 1;
    int i;


#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_create(%p, %d, %p, %p, %d, %p)", comm_old,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     ndims, dims, periods, reorder, comm_cart);
#endif
    if (!comm_old)
    {
        return tMPI_Error(comm_old, TMPI_ERR_COMM);
    }
    /* calculate the total number of procs in cartesian comm */
<<<<<<< HEAD
    for(i=0;i<ndims;i++)
=======
    for (i = 0; i < ndims; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        Ntot *= dims[i];
    }
    /* refuse to create if there's not enough procs */
    if (comm_old->grp.N < Ntot)
    {
<<<<<<< HEAD
        *comm_cart=TMPI_COMM_NULL;
=======
        *comm_cart = TMPI_COMM_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#if 1
        return tMPI_Error(comm_old, TMPI_ERR_CART_CREATE_NPROCS);
#endif
    }

    if (key >= Ntot)
<<<<<<< HEAD
        key=TMPI_UNDEFINED;
=======
    {
        key = TMPI_UNDEFINED;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (reorder)
    {
        tMPI_Cart_map(comm_old, ndims, dims, periods, &key);
    }

<<<<<<< HEAD
    if (key==TMPI_UNDEFINED)
    {
        color=TMPI_UNDEFINED;
=======
    if (key == TMPI_UNDEFINED)
    {
        color = TMPI_UNDEFINED;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    tMPI_Comm_split(comm_old, color, key, comm_cart);

    tMPI_Cart_init(comm_cart, ndims, dims, periods);

    return TMPI_SUCCESS;
}


int tMPI_Cart_sub(tMPI_Comm comm, int *remain_dims, tMPI_Comm *newcomm)
{
<<<<<<< HEAD
    int myrank;
    int ndims=0;
    int *dims=NULL;
    int *periods=NULL;
    int *oldcoords=NULL;
    int i;
    int ndims_notused=1;
    int color_notused=0;
=======
    int  myrank;
    int  ndims     = 0;
    int *dims      = NULL;
    int *periods   = NULL;
    int *oldcoords = NULL;
    int  i;
    int  ndims_notused = 1;
    int  color_notused = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Cart_sub(%p, %p, %p)", comm, remain_dims, newcomm);
#endif
    tMPI_Comm_rank(comm, &myrank);
<<<<<<< HEAD
    if ( comm->cart )
    {
        oldcoords=(int*)tMPI_Malloc(sizeof(int)*comm->cart->ndims);
        dims=(int*)tMPI_Malloc(sizeof(int)*comm->cart->ndims);
        periods=(int*)tMPI_Malloc(sizeof(int)*comm->cart->ndims);
=======
    if (comm->cart)
    {
        oldcoords = (int*)tMPI_Malloc(sizeof(int)*comm->cart->ndims);
        dims      = (int*)tMPI_Malloc(sizeof(int)*comm->cart->ndims);
        periods   = (int*)tMPI_Malloc(sizeof(int)*comm->cart->ndims);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        /* get old coordinates */
        tMPI_Cart_coords(comm, myrank, comm->cart->ndims, oldcoords);

<<<<<<< HEAD
        for(i=0;i<comm->cart->ndims;i++)
=======
        for (i = 0; i < comm->cart->ndims; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            if (remain_dims[i])
            {
                /* for the remaining dimensions, copy dimensionality data */
<<<<<<< HEAD
                dims[ndims]=comm->cart->dims[i];
                periods[ndims]=comm->cart->periods[i];
=======
                dims[ndims]    = comm->cart->dims[i];
                periods[ndims] = comm->cart->periods[i];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                ndims++;
            }
            else
            {
<<<<<<< HEAD
                /* base color on not used coordinates. We keep a 
=======
                /* base color on not used coordinates. We keep a
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                   ndims_notused index multiplier.*/
                color_notused += oldcoords[i]*ndims_notused;
                ndims_notused *= comm->cart->dims[i];
            }
        }
    }

    /* key=myrank, because we want the order to remain the same */
    tMPI_Comm_split(comm, color_notused, myrank, newcomm);
    tMPI_Cart_init(newcomm, ndims, dims, periods);

    if (oldcoords)
<<<<<<< HEAD
        free(oldcoords);
    if (dims)
        free(dims);
    if (periods)
        free(periods);

    return TMPI_SUCCESS;
}

=======
    {
        free(oldcoords);
    }
    if (dims)
    {
        free(dims);
    }
    if (periods)
    {
        free(periods);
    }

    return TMPI_SUCCESS;
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
