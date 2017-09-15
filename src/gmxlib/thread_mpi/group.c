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



/* Group query & manipulation functions */

tmpi_bool tMPI_In_group(tMPI_Group group)
{
<<<<<<< HEAD
    int i;
    struct tmpi_thread *cur;

    cur=tMPI_Get_current();
    for(i=0;i<group->N;i++)
    {
        if (group->peers[i] == cur)
            return TRUE;
=======
    int                 i;
    struct tmpi_thread *cur;

    cur = tMPI_Get_current();
    for (i = 0; i < group->N; i++)
    {
        if (group->peers[i] == cur)
        {
            return TRUE;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    return FALSE;
}

int tMPI_Group_size(tMPI_Group group, int *size)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Group_size(%p, %p)", group, size);
#endif

    if (group)
<<<<<<< HEAD
        *size = group->N;
    else
        *size = 0;
=======
    {
        *size = group->N;
    }
    else
    {
        *size = 0;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return TMPI_SUCCESS;
}

int tMPI_Group_rank(tMPI_Group group, int *rank)
<<<<<<< HEAD
{    
    int i;
=======
{
    int                 i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    struct tmpi_thread *cur;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Group_rank(%p, %p)", group, rank);
#endif
    if (!group)
<<<<<<< HEAD
        return TMPI_UNDEFINED;

    /* search for my id in the list of peers */
    cur=tMPI_Get_current();
    for(i=0;i<group->N;i++)
    {
        if (group->peers[i] == cur)
        {
            *rank=i;
=======
    {
        return TMPI_UNDEFINED;
    }

    /* search for my id in the list of peers */
    cur = tMPI_Get_current();
    for (i = 0; i < group->N; i++)
    {
        if (group->peers[i] == cur)
        {
            *rank = i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            return TMPI_SUCCESS;
        }
    }
    return TMPI_UNDEFINED;
}



tMPI_Group tMPI_Group_alloc(void)
{
    struct tmpi_group_ *ret;

<<<<<<< HEAD
    ret=(struct tmpi_group_*)tMPI_Malloc(sizeof(struct tmpi_group_));
    ret->peers=(struct tmpi_thread**)tMPI_Malloc(
                                sizeof(struct tmpi_thread*)*Nthreads);
    ret->N=0;
#if 0
    ret->Nrefs=1;
=======
    ret        = (struct tmpi_group_*)tMPI_Malloc(sizeof(struct tmpi_group_));
    ret->peers = (struct tmpi_thread**)tMPI_Malloc(
                sizeof(struct tmpi_thread*)*Nthreads);
    ret->N = 0;
#if 0
    ret->Nrefs = 1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

    return ret;
}

int tMPI_Group_free(tMPI_Group *group)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Group_free(%p)", group);
#endif
    if (group)
    {
        free((*group)->peers);
        free(*group);
    }
    return TMPI_SUCCESS;
}

int tMPI_Comm_group(tMPI_Comm comm, tMPI_Group *group)
{
<<<<<<< HEAD
    int i;
    struct tmpi_group_ *ret=tMPI_Group_alloc();
=======
    int                 i;
    struct tmpi_group_ *ret = tMPI_Group_alloc();
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_group(%p, %p)", comm, group);
#endif
<<<<<<< HEAD
    ret->N=comm->grp.N;
    for(i=0;i<comm->grp.N;i++)
    {
        ret->peers[i]=comm->grp.peers[i];
    }
    *group=ret;
#if 0
    if (comm)
    {
        *group=&(comm->grp);
    }
    else
    {
        *group=NULL;
=======
    ret->N = comm->grp.N;
    for (i = 0; i < comm->grp.N; i++)
    {
        ret->peers[i] = comm->grp.peers[i];
    }
    *group = ret;
#if 0
    if (comm)
    {
        *group = &(comm->grp);
    }
    else
    {
        *group = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
#endif

    return TMPI_SUCCESS;
}


int tMPI_Group_incl(tMPI_Group group, int n, int *ranks, tMPI_Group *newgroup)
{
<<<<<<< HEAD
    int i;
    tMPI_Group ng;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Group_incl(%p, %d, %p, %p, %p)", group, n, ranks, 
                       newgroup);
#endif
    /* just allocate and copy */
    ng=tMPI_Group_alloc();
    ng->N=n;
    for(i=0;i<n;i++)
=======
    int        i;
    tMPI_Group ng;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Group_incl(%p, %d, %p, %p, %p)", group, n, ranks,
                     newgroup);
#endif
    /* just allocate and copy */
    ng    = tMPI_Group_alloc();
    ng->N = n;
    for (i = 0; i < n; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (ranks[i] < 0 || !group || ranks[i] >= group->N)
        {
            return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_GROUP_RANK);
        }
<<<<<<< HEAD
        ng->peers[i]=group->peers[ranks[i]];
    }
    *newgroup=ng;
    return TMPI_SUCCESS;
}



=======
        ng->peers[i] = group->peers[ranks[i]];
    }
    *newgroup = ng;
    return TMPI_SUCCESS;
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
