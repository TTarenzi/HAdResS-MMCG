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


/*#define TMPI_DEBUG*/


/* If this is defined, thread_mpi will print a message when for every MPI
<<<<<<< HEAD
   call is called or returns. Useful for debugging MPI-related issues 
=======
   call is called or returns. Useful for debugging MPI-related issues
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   in the calling program. */
/*#define TMPI_TRACE*/

/* if this is defined, MPI will warn/hang/crash on practices that don't conform
   to the MPI standard (such as not calling tMPI_Comm_free on all threads that
   are part of the comm being freed). */
<<<<<<< HEAD
#define TMPI_STRICT 

/* whether to warn if there are mallocs at performance-critical sections 
=======
#define TMPI_STRICT

/* whether to warn if there are mallocs at performance-critical sections
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   (due to preallocations being too small) */
/*#define TMPI_WARN_MALLOC*/


/* the number of envelopes to allocate per thread-to-thread path */
#define N_EV_ALLOC 16

/* the normal maximum number of threads for pre-defined arrays
   (if the actual number of threads is bigger than this, it'll
    allocate/deallocate arrays, so no problems will arise). */
#define MAX_PREALLOC_THREADS 64

/* Whether to use lock-free lists using compare-and-swap (cmpxchg on x86)
<<<<<<< HEAD
   pointer functions. Message passing using blocking Send/Recv, and multicasts 
=======
   pointer functions. Message passing using blocking Send/Recv, and multicasts
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   are is still blocking, of course. */
#define TMPI_LOCK_FREE_LISTS

/* Whether to disable yielding to the OS scheduler during waits. Disabling
<<<<<<< HEAD
   this improves performance very slightly if Nthreads<=Ncores on an 
   otherwise idle system because waits have slightly lower latencies, but
   causes very poor performance if threads are competing for CPU time (for
   example, when Nthreads>Ncores, or another process is running on the 
   system.

   This option can be set with cmake. */
/*#define TMPI_WAIT_FOR_NO_ONE */ 



/* whether to enable double-copying (where the sender copies data to an 
   intermediate buffer for small enough buffers, allowing it to return
   from a blocking send call early. The receiver is free to copy from the
   original buffer while the sender is copying, possibly allowing them to
   work in parallel). 
   
=======
   this improves performance very slightly if Nthreads<=Ncores on an
   otherwise idle system because waits have slightly lower latencies, but
   causes very poor performance if threads are competing for CPU time (for
   example, when Nthreads>Ncores, or another process is running on the
   system.

   This option can be set with cmake. */
/*#define TMPI_WAIT_FOR_NO_ONE */



/* whether to enable double-copying (where the sender copies data to an
   intermediate buffer for small enough buffers, allowing it to return
   from a blocking send call early. The receiver is free to copy from the
   original buffer while the sender is copying, possibly allowing them to
   work in parallel).

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   This option can be set with cmake. */
/*#define TMPI_COPY_BUFFER*/


<<<<<<< HEAD
/* The size (in bytes) of the maximum transmission size for which double 
   copying is allowed (i.e. the sender doesn't wait for the receiver to 
=======
/* The size (in bytes) of the maximum transmission size for which double
   copying is allowed (i.e. the sender doesn't wait for the receiver to
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   become ready, but posts a copied buffer in its envelope).

   A size of 8192 bytes was chosen after some testing with Gromacs. */
#define COPY_BUFFER_SIZE 8192
#ifdef TMPI_COPY_BUFFER
/* We can separately specify whether we want copy buffers for send/recv or
   multicast communications: */
#define USE_SEND_RECV_COPY_BUFFER
#define USE_COLLECTIVE_COPY_BUFFER
#endif


/* The number of collective envelopes per comm object. This is the maximum
<<<<<<< HEAD
   number of simulataneous collective communications that can 
=======
   number of simulataneous collective communications that can
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   take place per comm object. If TMPI_NO_COPY_BUFFER is set, simultaneous
   collective communications don't happen and 2 is the right value.  */
#ifdef USE_COLLECTIVE_COPY_BUFFER
#define N_COLL_ENV 12
#else
#define N_COLL_ENV 2
#endif


<<<<<<< HEAD
/* Whether to do profiling of the number of MPI communication calls. A 
    report with the total number of calls for each communication function
    will be generated at MPI_Finalize(). 
    
=======
/* Whether to do profiling of the number of MPI communication calls. A
    report with the total number of calls for each communication function
    will be generated at MPI_Finalize().

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    This option can be set with cmake.*/
/*#define TMPI_PROFILE*/


<<<<<<< HEAD
/* whether to turn on thread affinity (required for NUMA optimizations) 
   if the number of threads to spawn is equal to the number of processors. */
#define TMPI_THREAD_AFFINITY


=======
/* whether to turn on thread affinity (required for NUMA optimizations)
   if the number of threads to spawn is equal to the number of processors. */
#define TMPI_THREAD_AFFINITY
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
