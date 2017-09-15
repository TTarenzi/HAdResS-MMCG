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
#include "p2p.h"


int tMPI_Wait(tMPI_Request *request, tMPI_Status *status)
{
<<<<<<< HEAD
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
    struct tmpi_req_ *rq;
=======
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
    struct tmpi_req_   *rq;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Wait(%p, %p)", request, status);
#endif
    if (!request || !(*request))
<<<<<<< HEAD
        return TMPI_SUCCESS;

    rq=*request;
    /* fix the pointers */
    rq->next=rq;
    rq->prev=rq;
=======
    {
        return TMPI_SUCCESS;
    }

    rq = *request;
    /* fix the pointers */
    rq->next = rq;
    rq->prev = rq;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* and wait for our request */
    do
    {
        if (tMPI_Test_single(cur, rq))
<<<<<<< HEAD
            break;
        tMPI_Wait_process_incoming(cur);
    } while(TRUE);

    rq->ev=NULL; /* we won't be using that envelope any more */
    ret=rq->error;
=======
        {
            break;
        }
        tMPI_Wait_process_incoming(cur);
    }
    while (TRUE);

    rq->ev = NULL; /* we won't be using that envelope any more */
    ret    = rq->error;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    tMPI_Set_status(rq, status);

    /* deallocate */
    tMPI_Return_req(rql, *request);

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Wait);
#endif
    return ret;
}

int tMPI_Test(tMPI_Request *request, int *flag, tMPI_Status *status)
{
<<<<<<< HEAD
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
    struct tmpi_req_ *rq;
=======
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
    struct tmpi_req_   *rq;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Test(%p, %p, %p)", request, flag, status);
#endif
    if (!request || !(*request))
<<<<<<< HEAD
        return TMPI_SUCCESS;

    rq=*request;
    /* fix the pointers */
    rq->next=rq;
    rq->prev=rq;

    /* and check our request */
    if (tMPI_Test_single(cur, rq))
        *flag=TRUE;

    ret=rq->error;
=======
    {
        return TMPI_SUCCESS;
    }

    rq = *request;
    /* fix the pointers */
    rq->next = rq;
    rq->prev = rq;

    /* and check our request */
    if (tMPI_Test_single(cur, rq))
    {
        *flag = TRUE;
    }

    ret = rq->error;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (rq->finished)
    {
        tMPI_Set_status(rq, status);
        /* deallocate */
        tMPI_Return_req(rql, *request);
<<<<<<< HEAD
        *request=TMPI_REQUEST_NULL;
=======
        *request = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Test);
#endif
    return ret;
}
















/* test multiple requests by first making a linked list of them, and then
   checking for their completion. Used in tMPI_{Test|Wait}{all|any|some}

<<<<<<< HEAD
   wait = whether to wait for incoming events 
   blocking = whether to block until all reqs are completed */
static void tMPI_Test_multi_req(struct tmpi_thread *cur, 
                                int count, tMPI_Request *array_of_requests,
                                tmpi_bool wait, tmpi_bool blocking)
{
    int i;
    struct tmpi_req_ *first=NULL, *last=NULL;

    /* construct the list of requests */
    for(i=0;i<count;i++)
    {
        struct tmpi_req_ *curr=array_of_requests[i];
        if (curr)
        {
            if (!first)
                first=curr;
=======
   wait = whether to wait for incoming events
   blocking = whether to block until all reqs are completed */
static void tMPI_Test_multi_req(struct tmpi_thread *cur,
                                int count, tMPI_Request *array_of_requests,
                                tmpi_bool wait, tmpi_bool blocking)
{
    int               i;
    struct tmpi_req_ *first = NULL, *last = NULL;

    /* construct the list of requests */
    for (i = 0; i < count; i++)
    {
        struct tmpi_req_ *curr = array_of_requests[i];
        if (curr)
        {
            if (!first)
            {
                first = curr;
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* fix the pointers */
            if (!last)
            {
                /* we connect to itself */
<<<<<<< HEAD
                last=curr;
                last->next=NULL;
                last->prev=NULL;
=======
                last       = curr;
                last->next = NULL;
                last->prev = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            else
            {
                /* we connect to the last */
<<<<<<< HEAD
                curr->next=NULL;
                curr->prev=last; 
                last->next=curr;
                last=curr;
=======
                curr->next = NULL;
                curr->prev = last;
                last->next = curr;
                last       = curr;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
    /* and wait for our request */
    do
    {
        if (tMPI_Test_multi(cur, first, NULL))
<<<<<<< HEAD
            break;
        if (wait)
            tMPI_Wait_process_incoming(cur);
    } while(blocking && wait);
=======
        {
            break;
        }
        if (wait)
        {
            tMPI_Wait_process_incoming(cur);
        }
    }
    while (blocking && wait);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}



int tMPI_Waitall(int count, tMPI_Request *array_of_requests,
                 tMPI_Status *array_of_statuses)
{
<<<<<<< HEAD
    int i;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
=======
    int                 i;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
<<<<<<< HEAD
    tMPI_Trace_print("tMPI_Waitall(%d, %p, %p)", count, array_of_requests, 
                       array_of_statuses);
=======
    tMPI_Trace_print("tMPI_Waitall(%d, %p, %p)", count, array_of_requests,
                     array_of_statuses);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    tMPI_Test_multi_req(cur, count, array_of_requests, TRUE, TRUE);

    /* deallocate the now finished requests */
<<<<<<< HEAD
    for(i=0;i<count;i++)
=======
    for (i = 0; i < count; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (array_of_requests[i])
        {
            if (array_of_statuses)
            {
                tMPI_Set_status(array_of_requests[i], &(array_of_statuses[i]));
            }
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
<<<<<<< HEAD
                ret=TMPI_ERR_IN_STATUS; 
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i]=TMPI_REQUEST_NULL;
=======
                ret = TMPI_ERR_IN_STATUS;
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i] = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Waitall);
#endif
    return ret;
}

int tMPI_Testall(int count, tMPI_Request *array_of_requests,
                 int *flag, tMPI_Status *array_of_statuses)
{
<<<<<<< HEAD
    int i;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
=======
    int                 i;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
<<<<<<< HEAD
    tMPI_Trace_print("tMPI_Testall(%d, %p, %p, %p)", count, array_of_requests, 
                       flag, array_of_statuses);
=======
    tMPI_Trace_print("tMPI_Testall(%d, %p, %p, %p)", count, array_of_requests,
                     flag, array_of_statuses);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    tMPI_Test_multi_req(cur, count, array_of_requests, FALSE, TRUE);

    if (flag)
<<<<<<< HEAD
        *flag=1;
    /* deallocate the possibly finished requests */
    for(i=0;i<count;i++)
=======
    {
        *flag = 1;
    }
    /* deallocate the possibly finished requests */
    for (i = 0; i < count; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (array_of_requests[i] && array_of_requests[i]->finished)

        {
            if (array_of_statuses)
            {
                tMPI_Set_status(array_of_requests[i], &(array_of_statuses[i]));
            }
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
<<<<<<< HEAD
                ret=TMPI_ERR_IN_STATUS; 
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i]=TMPI_REQUEST_NULL;
=======
                ret = TMPI_ERR_IN_STATUS;
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i] = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        else
        {
            if (flag)
<<<<<<< HEAD
                *flag=0;
        } 
=======
            {
                *flag = 0;
            }
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Testall);
#endif
    return ret;
}


int tMPI_Waitany(int count, tMPI_Request *array_of_requests, int *index,
                 tMPI_Status *status)
{
<<<<<<< HEAD
    int i;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
=======
    int                 i;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
<<<<<<< HEAD
    tMPI_Trace_print("tMPI_Waitany(%d, %p, %p, %p)", count, array_of_requests, 
                       index, status);
=======
    tMPI_Trace_print("tMPI_Waitany(%d, %p, %p, %p)", count, array_of_requests,
                     index, status);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

    tMPI_Test_multi_req(cur, count, array_of_requests, TRUE, FALSE);

    /* deallocate the possibly finished requests */
<<<<<<< HEAD
    for(i=0;i<count;i++)
=======
    for (i = 0; i < count; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (array_of_requests[i] && array_of_requests[i]->finished)
        {
            tMPI_Set_status(array_of_requests[i], status);
            if (index)
<<<<<<< HEAD
                *index=i;
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
                ret=TMPI_ERR_IN_STATUS; 
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i]=TMPI_REQUEST_NULL;
=======
            {
                *index = i;
            }
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
                ret = TMPI_ERR_IN_STATUS;
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i] = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* and we only need one */
            break;
        }
    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Waitany);
#endif
    return ret;
}



int tMPI_Testany(int count, tMPI_Request *array_of_requests, int *index,
                 int *flag, tMPI_Status *status)
{
<<<<<<< HEAD
    int i;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
=======
    int                 i;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
<<<<<<< HEAD
    tMPI_Trace_print("tMPI_Testany(%d, %p, %p %p, %p)", count, 
=======
    tMPI_Trace_print("tMPI_Testany(%d, %p, %p %p, %p)", count,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     array_of_requests, flag, index, status);
#endif

    tMPI_Test_multi_req(cur, count, array_of_requests, FALSE, FALSE);

    if (flag)
<<<<<<< HEAD
        *flag=0;
    if (index)
        *index=TMPI_UNDEFINED;
    /* deallocate the possibly finished requests */
    for(i=0;i<count;i++)
=======
    {
        *flag = 0;
    }
    if (index)
    {
        *index = TMPI_UNDEFINED;
    }
    /* deallocate the possibly finished requests */
    for (i = 0; i < count; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (array_of_requests[i] && array_of_requests[i]->finished)
        {
            tMPI_Set_status(array_of_requests[i], status);
            if (index)
<<<<<<< HEAD
                *index=i;
            if (flag)
                *flag=1;
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
                ret=TMPI_ERR_IN_STATUS; 
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i]=TMPI_REQUEST_NULL;
=======
            {
                *index = i;
            }
            if (flag)
            {
                *flag = 1;
            }
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
                ret = TMPI_ERR_IN_STATUS;
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i] = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            /* and we only need one */
            break;
        }
    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Testany);
#endif
    return ret;
}



int tMPI_Waitsome(int incount, tMPI_Request *array_of_requests,
                  int *outcount, int *array_of_indices,
                  tMPI_Status *array_of_statuses)
{
<<<<<<< HEAD
    int i;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
=======
    int                 i;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
<<<<<<< HEAD
    tMPI_Trace_print("tMPI_Waitsome(%d, %p, %p, %p, %p)", incount, 
                     array_of_requests, outcount, array_of_indices, 
=======
    tMPI_Trace_print("tMPI_Waitsome(%d, %p, %p, %p, %p)", incount,
                     array_of_requests, outcount, array_of_indices,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     array_of_statuses);
#endif
    tMPI_Test_multi_req(cur, incount, array_of_requests, TRUE, FALSE);

<<<<<<< HEAD
    (*outcount)=0;
    /* deallocate the possibly finished requests */
    for(i=0;i<incount;i++)
=======
    (*outcount) = 0;
    /* deallocate the possibly finished requests */
    for (i = 0; i < incount; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (array_of_requests[i] && array_of_requests[i]->finished)
        {
            array_of_indices[*outcount]++;
            (*outcount)++;
            if (array_of_statuses)
            {
                tMPI_Set_status(array_of_requests[i], &(array_of_statuses[i]));
            }
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
<<<<<<< HEAD
                ret=TMPI_ERR_IN_STATUS; 
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i]=TMPI_REQUEST_NULL;
=======
                ret = TMPI_ERR_IN_STATUS;
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i] = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Waitsome);
#endif
    return ret;
}

int tMPI_Testsome(int incount, tMPI_Request *array_of_requests,
                  int *outcount, int *array_of_indices,
                  tMPI_Status *array_of_statuses)
{
<<<<<<< HEAD
    int i;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();
    struct req_list *rql=&(cur->rql);
=======
    int                 i;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
<<<<<<< HEAD
    tMPI_Trace_print("tMPI_Testsome(%d, %p, %p, %p, %p)", incount, 
                     array_of_requests, outcount, array_of_indices, 
=======
    tMPI_Trace_print("tMPI_Testsome(%d, %p, %p, %p, %p)", incount,
                     array_of_requests, outcount, array_of_indices,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     array_of_statuses);
#endif
    tMPI_Test_multi_req(cur, incount, array_of_requests, FALSE, TRUE);

<<<<<<< HEAD
    (*outcount)=0;
    /* deallocate the possibly finished requests */
    for(i=0;i<incount;i++)
=======
    (*outcount) = 0;
    /* deallocate the possibly finished requests */
    for (i = 0; i < incount; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (array_of_requests[i] && array_of_requests[i]->finished)
        {
            array_of_indices[*outcount]++;
            (*outcount)++;
            if (array_of_statuses)
            {
                tMPI_Set_status(array_of_requests[i], &(array_of_statuses[i]));
            }
            if (array_of_requests[i]->error != TMPI_SUCCESS)
            {
<<<<<<< HEAD
                ret=TMPI_ERR_IN_STATUS; 
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i]=TMPI_REQUEST_NULL;
=======
                ret = TMPI_ERR_IN_STATUS;
            }
            tMPI_Return_req(rql, array_of_requests[i]);
            array_of_requests[i] = TMPI_REQUEST_NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }


#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Testsome);
#endif
    return ret;
}
<<<<<<< HEAD


=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
