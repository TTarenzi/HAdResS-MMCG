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



/* Include the defines that determine which thread library to use.
 * We do not use HAVE_PTHREAD_H directly, since we might want to
 * turn off thread support explicity (e.g. for debugging).
 */

#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef THREAD_WINDOWS

/* the win32 header */
#include <windows.h>


#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


#include "thread_mpi/atomic.h"
#include "thread_mpi/threads.h"
#include "impl.h"

#include "winthreads.h"

<<<<<<< HEAD
/*! \brief System mutex for all one-time initialization 
 *
 *  This static variable is necessary in order to make the header file 
=======
/*! \brief System mutex for all one-time initialization
 *
 *  This static variable is necessary in order to make the header file
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 */

static CRITICAL_SECTION mutex_init;   /* mutex for initializing mutexes */
static CRITICAL_SECTION once_init;    /* mutex for initializing barriers */
static CRITICAL_SECTION cond_init;    /* mutex for initializing thread_conds */
static CRITICAL_SECTION barrier_init; /* mutex for initializing barriers */


/* spinlock for initializing the above mutexes */
<<<<<<< HEAD
static tMPI_Spinlock_t init_init=TMPI_SPINLOCK_INITIALIZER;

/* whether tMPI_Thread_create has initialized these mutexes */
static tMPI_Atomic_t init_inited={ 0 };

/* whether the main thread affinity has been set */
static tMPI_Spinlock_t main_thread_aff_lock=TMPI_SPINLOCK_INITIALIZER;
static tMPI_Atomic_t main_thread_aff_set={ 0 };
=======
static tMPI_Spinlock_t init_init = TMPI_SPINLOCK_INITIALIZER;

/* whether tMPI_Thread_create has initialized these mutexes */
static tMPI_Atomic_t init_inited = { 0 };

/* whether the main thread affinity has been set */
static tMPI_Spinlock_t main_thread_aff_lock = TMPI_SPINLOCK_INITIALIZER;
static tMPI_Atomic_t   main_thread_aff_set  = { 0 };

/* mutex for managing  thread IDs */
static CRITICAL_SECTION thread_id_list_lock;
typedef struct
{
    DWORD               thread_id; /* the thread ID as returned by GetCurrentTreadID() */
    struct tMPI_Thread* th;        /* the associated tMPI thread structure */
} thread_id_list_t;
/* the size of the thrread id list */
static int               Nalloc_thread_id_list = 0;
/* the number of elements in the thread id list */
static int               N_thread_id_list = 0;
/* the thread ID list */
static thread_id_list_t *thread_id_list;



/* data structure to keep track of thread key destructors. */
typedef struct
{
    void (*destructor)(void*);
    DWORD key;
} thread_key_destructors;

static thread_key_destructors *destructors = NULL;


>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/*
    NUMA and Processor Group awareness support.

<<<<<<< HEAD
    NUMA support is implemented to maximize the chance that memory access 
    patterns remain Local to the NUMA node.
    NUMA node processor affinity is utilized to prevent scheduler associated 
    drift across NUMA nodes.
    Process Group support is implemented to enable > 64 processors to be 
    utilized.  This is only supported when building 64bit.

    The high level approach is:
    1. Build a description of CPU topology, including processor numbers, NUMA 
        node numbers, and affinity masks.
    2. For processor intensive worker threads, create threads such that 
        the processor affinity and thread stack is kept local within a NUMA node.
    3. Employ simple round-robin affinity and node assignment approach when 
        creating threads.
    4. Use GetProcAddress() to obtain function pointers to functions that 
        are operating system version dependent, to allow maximum binary 
        compatibility. 

    Scott Field (sfield@microsoft.com)      Jan-2011    
*/
=======
    NUMA support is implemented to maximize the chance that memory access
    patterns remain Local to the NUMA node.
    NUMA node processor affinity is utilized to prevent scheduler associated
    drift across NUMA nodes.
    Process Group support is implemented to enable > 64 processors to be
    utilized.  This is only supported when building 64bit.

    The high level approach is:
    1. Build a description of CPU topology, including processor numbers, NUMA
        node numbers, and affinity masks.
    2. For processor intensive worker threads, create threads such that
        the processor affinity and thread stack is kept local within a NUMA node.
    3. Employ simple round-robin affinity and node assignment approach when
        creating threads.
    4. Use GetProcAddress() to obtain function pointers to functions that
        are operating system version dependent, to allow maximum binary
        compatibility.

    Scott Field (sfield@microsoft.com)      Jan-2011
 */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


typedef struct {
    PROCESSOR_NUMBER    ProcessorNumber;
    GROUP_AFFINITY      GroupAffinity;
    USHORT              NumaNodeNumber;
} MPI_NUMA_PROCESSOR_INFO;


/* thread/processor index, to allow setting round-robin affinity. */
<<<<<<< HEAD
volatile ULONG g_ulThreadIndex;                 
/* a value of zero implies the system is not NUMA */
ULONG g_ulHighestNumaNodeNumber=0;
/* total number of processors in g_MPI_ProcessInfo array */
ULONG g_ulTotalProcessors;
/* array describing available processors, affinity masks, and NUMA node */
MPI_NUMA_PROCESSOR_INFO *g_MPI_ProcessorInfo;   

/* function prototypes and variables to support obtaining function addresses 
   dynamically -- supports down-level operating systems */

typedef BOOL (WINAPI *func_GetNumaHighestNodeNumber_t)( PULONG 
                                                        HighestNodeNumber );
typedef BOOL (WINAPI *func_SetThreadGroupAffinity_t)( HANDLE hThread, 
                            const GROUP_AFFINITY *GroupAffinity, 
                            PGROUP_AFFINITY PreviousGroupAffinity );
typedef BOOL (WINAPI *func_SetThreadIdealProcessorEx_t)( HANDLE hThread, 
                            PPROCESSOR_NUMBER lpIdealProcessor, 
                            PPROCESSOR_NUMBER lpPreviousIdealProcessor );
typedef BOOL (WINAPI *func_GetNumaNodeProcessorMaskEx_t)( USHORT Node, 
                            PGROUP_AFFINITY ProcessorMask );
typedef BOOL (WINAPI *func_GetNumaProcessorNodeEx_t)( 
                            PPROCESSOR_NUMBER Processor, 
                            PUSHORT NodeNumber );
typedef VOID (WINAPI *func_GetCurrentProcessorNumberEx_t)( 
                            PPROCESSOR_NUMBER ProcNumber );

typedef HANDLE (WINAPI *func_CreateRemoteThreadEx_t)(
                            HANDLE hProcess,
                            LPSECURITY_ATTRIBUTES lpThreadAttributes,
                            SIZE_T dwStackSize,
                            LPTHREAD_START_ROUTINE lpStartAddress,
                            LPVOID lpParameter,
                            DWORD dwCreationFlags,
                            LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList,
                            LPDWORD lpThreadId);

typedef BOOL (WINAPI *func_InitializeProcThreadAttributeList_t)(
                            LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList, 
                            DWORD dwAttributeCount, 
                            DWORD dwFlags, 
                            PSIZE_T lpSize);
typedef BOOL (WINAPI *func_UpdateProcThreadAttribute_t)(
                            LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList,
                            DWORD dwFlags,
                            DWORD_PTR Attribute,
                            PVOID lpValue,
                            SIZE_T cbSize,
                            PVOID lpPreviousValue,
                            PSIZE_T lpReturnSize);
typedef VOID (WINAPI *func_DeleteProcThreadAttributeList_t)(
                            LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList);
=======
volatile ULONG           g_ulThreadIndex;
/* a value of zero implies the system is not NUMA */
ULONG                    g_ulHighestNumaNodeNumber = 0;
/* total number of processors in g_MPI_ProcessInfo array */
ULONG                    g_ulTotalProcessors;
/* array describing available processors, affinity masks, and NUMA node */
MPI_NUMA_PROCESSOR_INFO *g_MPI_ProcessorInfo = NULL;

/* function prototypes and variables to support obtaining function addresses
   dynamically -- supports down-level operating systems */

typedef BOOL (WINAPI *func_GetNumaHighestNodeNumber_t)( PULONG
                                                        HighestNodeNumber );
typedef DWORD (WINAPI *func_SetThreadIdealProcessor_t)( HANDLE hThread,
                                                        DWORD dwIdealProcessor );
typedef BOOL (WINAPI *func_SetThreadGroupAffinity_t)( HANDLE hThread,
                                                      const GROUP_AFFINITY *GroupAffinity,
                                                      PGROUP_AFFINITY PreviousGroupAffinity );
typedef BOOL (WINAPI *func_SetThreadIdealProcessorEx_t)( HANDLE hThread,
                                                         PPROCESSOR_NUMBER lpIdealProcessor,
                                                         PPROCESSOR_NUMBER lpPreviousIdealProcessor );
typedef BOOL (WINAPI *func_GetNumaNodeProcessorMaskEx_t)( USHORT Node,
                                                          PGROUP_AFFINITY ProcessorMask );
typedef BOOL (WINAPI *func_GetNumaProcessorNodeEx_t)(
        PPROCESSOR_NUMBER Processor,
        PUSHORT NodeNumber );
typedef VOID (WINAPI *func_GetCurrentProcessorNumberEx_t)(
        PPROCESSOR_NUMBER ProcNumber );

typedef HANDLE (WINAPI *func_CreateRemoteThreadEx_t)(
        HANDLE hProcess,
        LPSECURITY_ATTRIBUTES lpThreadAttributes,
        SIZE_T dwStackSize,
        LPTHREAD_START_ROUTINE lpStartAddress,
        LPVOID lpParameter,
        DWORD dwCreationFlags,
        LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList,
        LPDWORD lpThreadId);

typedef BOOL (WINAPI *func_InitializeProcThreadAttributeList_t)(
        LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList,
        DWORD dwAttributeCount,
        DWORD dwFlags,
        PSIZE_T lpSize);
typedef BOOL (WINAPI *func_UpdateProcThreadAttribute_t)(
        LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList,
        DWORD dwFlags,
        DWORD_PTR Attribute,
        PVOID lpValue,
        SIZE_T cbSize,
        PVOID lpPreviousValue,
        PSIZE_T lpReturnSize);
typedef VOID (WINAPI *func_DeleteProcThreadAttributeList_t)(
        LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
typedef DWORD (WINAPI *func_GetActiveProcessorCount_t)(WORD GroupNumber);
typedef WORD (WINAPI *func_GetActiveProcessorGroupCount_t)(void);


/* WinXP SP2, WinXP64, WinSrv 2003 */
<<<<<<< HEAD
func_GetNumaHighestNodeNumber_t             func_GetNumaHighestNodeNumber;              
/* Windows 7, WinSrv 2008R2 */
func_SetThreadGroupAffinity_t               func_SetThreadGroupAffinity;                
func_SetThreadIdealProcessorEx_t            func_SetThreadIdealProcessorEx; 
func_GetNumaNodeProcessorMaskEx_t           func_GetNumaNodeProcessorMaskEx;
func_GetNumaProcessorNodeEx_t               func_GetNumaProcessorNodeEx;    
func_GetCurrentProcessorNumberEx_t          func_GetCurrentProcessorNumberEx;
func_GetActiveProcessorCount_t              func_GetActiveProcessorCount;    
func_GetActiveProcessorGroupCount_t         func_GetActiveProcessorGroupCount;
func_CreateRemoteThreadEx_t                 func_CreateRemoteThreadEx;        
/* Windows Vista, WinSrv 2008 */ 
func_InitializeProcThreadAttributeList_t    func_InitializeProcThreadAttributeList;     
func_UpdateProcThreadAttribute_t            func_UpdateProcThreadAttribute;             
func_DeleteProcThreadAttributeList_t        func_DeleteProcThreadAttributeList;         


/* Set the main thread's affinity */
static int tMPI_Set_main_thread_affinity(void)
{
    /* calling thread PROCESSOR_NUMBER */
    PROCESSOR_NUMBER CurrentProcessorNumber;      
    /* calling thread GROUP_AFFINITY */
    GROUP_AFFINITY CurrentThreadGroupAffinity; 
    /* calling thread NUMA node */
    USHORT CurrentNumaNodeNumber;


    /* we can pre-check because it's atomic */
    if (tMPI_Atomic_get(&main_thread_aff_set) == 0)
    {
        /* this can be a spinlock because the chances of collision are low. */
        tMPI_Spinlock_lock( &main_thread_aff_lock );
        if( g_ulHighestNumaNodeNumber != 0 )
        {
            func_GetCurrentProcessorNumberEx(&CurrentProcessorNumber);


            /* set the NUMA node affinity for the current thread
               failures to set the current thread affinity are ignored, 
               as a fringe case can arise on >32 processor systems with a 32bit 
               build/code.
               */
            func_SetThreadIdealProcessorEx(GetCurrentThread(), 
                                           &CurrentProcessorNumber, 
                                           NULL);

            if(func_GetNumaProcessorNodeEx(&CurrentProcessorNumber, 
                                           &CurrentNumaNodeNumber))
            {
                /* for the NUMA node number associated with the current processor 
                   number, get the group affinity mask */
                if(func_GetNumaNodeProcessorMaskEx(CurrentNumaNodeNumber, 
                                                   &CurrentThreadGroupAffinity))
                {
                    /* set the current thread affinity to prevent it from running on 
                       other NUMA nodes */
                    func_SetThreadGroupAffinity(GetCurrentThread(), 
                                                &CurrentThreadGroupAffinity, 
                                                NULL);
                }
            }
        }
        else
        {
            /* No NUMA. For now, we just do a similar thing. */
            if ( (func_GetCurrentProcessorNumberEx != NULL)  &&
                 (func_SetThreadIdealProcessorEx))
            {
                func_GetCurrentProcessorNumberEx(&CurrentProcessorNumber);
                func_SetThreadIdealProcessorEx(GetCurrentThread(), 
                                               &CurrentProcessorNumber, 
                                               NULL);
            }
        }
        tMPI_Atomic_set( &main_thread_aff_set, 1);
        tMPI_Spinlock_unlock( &main_thread_aff_lock );
    }
    return 0;
}

/*  returns 0 on success.
    Success is returned if the system is non-NUMA, OR the system doesn't 
    support appropriate NUMA APIs, OR the system is NUMA and we successfully 
    initialized support.
    
    returns -1 on error.
    This can happen if an API returned an error, a memory allocation failed, or 
    we failed to initialize affinity mapping information.
*/
int tMPI_Init_NUMA(void)
{
    /* module handle to kernel32.dll -- we already reference it, so it's already loaded */
    HMODULE hModKernel32 = NULL;                    
    /* 0-based NUMA node count -- does not imply all nodes have available (eg: hot-plug) processors */
    ULONG ulHighestNumaNodeNumber;                  
    /* total number of processors available per affinity masks */
    DWORD dwTotalProcessors = 0;                    
    ULONG i = 0;

    /* calling thread PROCESSOR_NUMBER */
    PROCESSOR_NUMBER CurrentProcessorNumber;      
    /* calling thread GROUP_AFFINITY */
    GROUP_AFFINITY CurrentThreadGroupAffinity; 
    /* calling thread NUMA node */
    USHORT CurrentNumaNodeNumber;
=======
func_GetNumaHighestNodeNumber_t             func_GetNumaHighestNodeNumber;
func_SetThreadIdealProcessor_t              func_SetThreadIdealProcessor;
/* Windows 7, WinSrv 2008R2 */
func_SetThreadGroupAffinity_t               func_SetThreadGroupAffinity;
func_SetThreadIdealProcessorEx_t            func_SetThreadIdealProcessorEx;
func_GetNumaNodeProcessorMaskEx_t           func_GetNumaNodeProcessorMaskEx;
func_GetNumaProcessorNodeEx_t               func_GetNumaProcessorNodeEx;
func_GetCurrentProcessorNumberEx_t          func_GetCurrentProcessorNumberEx;
func_GetActiveProcessorCount_t              func_GetActiveProcessorCount;
func_GetActiveProcessorGroupCount_t         func_GetActiveProcessorGroupCount;
func_CreateRemoteThreadEx_t                 func_CreateRemoteThreadEx;
/* Windows Vista, WinSrv 2008 */
func_InitializeProcThreadAttributeList_t    func_InitializeProcThreadAttributeList;
func_UpdateProcThreadAttribute_t            func_UpdateProcThreadAttribute;
func_DeleteProcThreadAttributeList_t        func_DeleteProcThreadAttributeList;



/*  returns 0 on success.
    Success is returned if the system is non-NUMA, OR the system doesn't
    support appropriate NUMA APIs, OR the system is NUMA and we successfully
    initialized support.

    returns -1 on error.
    This can happen if an API returned an error, a memory allocation failed, or
    we failed to initialize affinity mapping information.
 */
int tMPI_Init_NUMA(void)
{
    /* module handle to kernel32.dll -- we already reference it, so it's already loaded */
    HMODULE hModKernel32 = NULL;
    /* 0-based NUMA node count -- does not imply all nodes have available (eg: hot-plug) processors */
    ULONG   ulHighestNumaNodeNumber;
    /* total number of processors available per affinity masks */
    DWORD   dwTotalProcessors = 0;
    ULONG   i                 = 0;

    /* calling thread PROCESSOR_NUMBER */
    PROCESSOR_NUMBER CurrentProcessorNumber;
    /* calling thread GROUP_AFFINITY */
    /*GROUP_AFFINITY CurrentThreadGroupAffinity; */
    /* calling thread NUMA node */
    /*USHORT CurrentNumaNodeNumber;*/
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    WORD wActiveGroupCount;
    WORD GroupIndex;

    /* array of processor information structures */
<<<<<<< HEAD
    MPI_NUMA_PROCESSOR_INFO *pMPI_ProcessorInfo = NULL; 
=======
    MPI_NUMA_PROCESSOR_INFO *pMPI_ProcessorInfo = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* assume an error condition */
    int iRet = -1;

    hModKernel32 = GetModuleHandleA("kernel32.dll");

<<<<<<< HEAD
    if( hModKernel32 == NULL )
=======
    if (hModKernel32 == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return 0;
    }

<<<<<<< HEAD
    /* obtain addresses of relevant NUMA functions, most of which are 
       Windows 7 / Windows Server 2008R2 only functions
       this is done using GetProcAddress to enable the binary to run on older 
       Windows versions.
    */

    func_GetNumaHighestNodeNumber = (func_GetNumaHighestNodeNumber_t) GetProcAddress( hModKernel32, "GetNumaHighestNodeNumber" );

    if( func_GetNumaHighestNodeNumber == NULL )
=======
    /* obtain addresses of relevant NUMA functions, most of which are
       Windows 7 / Windows Server 2008R2 only functions
       this is done using GetProcAddress to enable the binary to run on older
       Windows versions.
     */

    func_GetNumaHighestNodeNumber = (func_GetNumaHighestNodeNumber_t) GetProcAddress( hModKernel32, "GetNumaHighestNodeNumber" );
    func_SetThreadIdealProcessor  = (func_SetThreadIdealProcessor_t) GetProcAddress( hModKernel32, "SetThreadIdealProcessor" );

    if (func_GetNumaHighestNodeNumber == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return 0;
    }

<<<<<<< HEAD
    /* determine if we're on a NUMA system and if so, determine the number of 
       (potential) nodes */

    if(!func_GetNumaHighestNodeNumber( &ulHighestNumaNodeNumber ))
=======
    /* determine if we're on a NUMA system and if so, determine the number of
       (potential) nodes */

    if (!func_GetNumaHighestNodeNumber( &ulHighestNumaNodeNumber ))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return -1;
    }

<<<<<<< HEAD
    if( ulHighestNumaNodeNumber == 0 )
    {
        /* system is not NUMA */
        return 0;
    }

    func_SetThreadGroupAffinity = (func_SetThreadGroupAffinity_t)GetProcAddress( hModKernel32, "SetThreadGroupAffinity" );
    func_SetThreadIdealProcessorEx = (func_SetThreadIdealProcessorEx_t)GetProcAddress( hModKernel32, "SetThreadIdealProcessorEx" );
    func_CreateRemoteThreadEx = (func_CreateRemoteThreadEx_t)GetProcAddress( hModKernel32, "CreateRemoteThreadEx" );
    func_GetNumaNodeProcessorMaskEx = (func_GetNumaNodeProcessorMaskEx_t)GetProcAddress( hModKernel32, "GetNumaNodeProcessorMaskEx" );
    func_GetNumaProcessorNodeEx = (func_GetNumaProcessorNodeEx_t)GetProcAddress( hModKernel32, "GetNumaProcessorNodeEx" );
    func_GetCurrentProcessorNumberEx = (func_GetCurrentProcessorNumberEx_t)GetProcAddress( hModKernel32, "GetCurrentProcessorNumberEx" );
    func_GetActiveProcessorCount = (func_GetActiveProcessorCount_t)GetProcAddress( hModKernel32, "GetActiveProcessorCount" );
    func_GetActiveProcessorGroupCount = (func_GetActiveProcessorGroupCount_t)GetProcAddress( hModKernel32, "GetActiveProcessorGroupCount" );
    func_InitializeProcThreadAttributeList = (func_InitializeProcThreadAttributeList_t)GetProcAddress( hModKernel32, "InitializeProcThreadAttributeList" );
    func_UpdateProcThreadAttribute = (func_UpdateProcThreadAttribute_t)GetProcAddress( hModKernel32, "UpdateProcThreadAttribute" );
    func_DeleteProcThreadAttributeList = (func_DeleteProcThreadAttributeList_t)GetProcAddress( hModKernel32, "DeleteProcThreadAttributeList" );

    if( (func_SetThreadGroupAffinity == NULL) ||
        (func_SetThreadIdealProcessorEx == NULL) ||
        (func_CreateRemoteThreadEx == NULL) ||
        (func_GetNumaNodeProcessorMaskEx == NULL) ||
        (func_GetNumaProcessorNodeEx == NULL) ||
        (func_GetCurrentProcessorNumberEx == NULL) ||
        (func_GetActiveProcessorCount == NULL) ||
        (func_GetActiveProcessorGroupCount == NULL) ||
        (func_InitializeProcThreadAttributeList == NULL) ||
        (func_UpdateProcThreadAttribute == NULL) ||
        (func_DeleteProcThreadAttributeList == NULL) )
    {
        /* if any addresses couldn't be located, assume NUMA functionality 
           isn't supported */
        return 0;
    }
=======


    func_SetThreadGroupAffinity            = (func_SetThreadGroupAffinity_t)GetProcAddress( hModKernel32, "SetThreadGroupAffinity" );
    func_SetThreadIdealProcessorEx         = (func_SetThreadIdealProcessorEx_t)GetProcAddress( hModKernel32, "SetThreadIdealProcessorEx" );
    func_CreateRemoteThreadEx              = (func_CreateRemoteThreadEx_t)GetProcAddress( hModKernel32, "CreateRemoteThreadEx" );
    func_GetNumaNodeProcessorMaskEx        = (func_GetNumaNodeProcessorMaskEx_t)GetProcAddress( hModKernel32, "GetNumaNodeProcessorMaskEx" );
    func_GetNumaProcessorNodeEx            = (func_GetNumaProcessorNodeEx_t)GetProcAddress( hModKernel32, "GetNumaProcessorNodeEx" );
    func_GetCurrentProcessorNumberEx       = (func_GetCurrentProcessorNumberEx_t)GetProcAddress( hModKernel32, "GetCurrentProcessorNumberEx" );
    func_GetActiveProcessorCount           = (func_GetActiveProcessorCount_t)GetProcAddress( hModKernel32, "GetActiveProcessorCount" );
    func_GetActiveProcessorGroupCount      = (func_GetActiveProcessorGroupCount_t)GetProcAddress( hModKernel32, "GetActiveProcessorGroupCount" );
    func_InitializeProcThreadAttributeList = (func_InitializeProcThreadAttributeList_t)GetProcAddress( hModKernel32, "InitializeProcThreadAttributeList" );
    func_UpdateProcThreadAttribute         = (func_UpdateProcThreadAttribute_t)GetProcAddress( hModKernel32, "UpdateProcThreadAttribute" );
    func_DeleteProcThreadAttributeList     = (func_DeleteProcThreadAttributeList_t)GetProcAddress( hModKernel32, "DeleteProcThreadAttributeList" );

    if ( (func_SetThreadGroupAffinity == NULL) ||
         (func_SetThreadIdealProcessorEx == NULL) ||
         (func_CreateRemoteThreadEx == NULL) ||
         (func_GetNumaNodeProcessorMaskEx == NULL) ||
         (func_GetNumaProcessorNodeEx == NULL) ||
         (func_GetCurrentProcessorNumberEx == NULL) ||
         (func_GetActiveProcessorCount == NULL) ||
         (func_GetActiveProcessorGroupCount == NULL) ||
         (func_InitializeProcThreadAttributeList == NULL) ||
         (func_UpdateProcThreadAttribute == NULL) ||
         (func_DeleteProcThreadAttributeList == NULL) )
    {
        /* if any addresses couldn't be located, assume NUMA functionality
           isn't supported */
        return 0;
    }
#if 0
    if (ulHighestNumaNodeNumber == 0)
    {
        /* system is not NUMA */
        return 0;
    }
#endif
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* count the active processors across the groups */

    func_GetCurrentProcessorNumberEx(&CurrentProcessorNumber);

    wActiveGroupCount = func_GetActiveProcessorGroupCount();
<<<<<<< HEAD
    
    dwTotalProcessors = func_GetActiveProcessorCount( ALL_PROCESSOR_GROUPS );

#if !((defined WIN64 || defined _WIN64))
    /* WOW64 doesn't allow setting the affinity correctly beyond 32 
       processors -- the KAFFINITY mask is only 32 bits wide
       This check is only here for completeness -- large systems should be 
       running 64bit Gromacs code, where the processor quantity is not 
       constrained.
       By failing here, the WOW64 32bit client will use normal CreateThread(), 
       which can schedule up to 64 un-affinitized threads
    */

    if( dwTotalProcessors > 32 )
=======

    dwTotalProcessors = func_GetActiveProcessorCount( ALL_PROCESSOR_GROUPS );

#if !((defined WIN64 || defined _WIN64))
    /* WOW64 doesn't allow setting the affinity correctly beyond 32
       processors -- the KAFFINITY mask is only 32 bits wide
       This check is only here for completeness -- large systems should be
       running 64bit Gromacs code, where the processor quantity is not
       constrained.
       By failing here, the WOW64 32bit client will use normal CreateThread(),
       which can schedule up to 64 un-affinitized threads
     */

    if (dwTotalProcessors > 32)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return 0;
    }
#endif

    /* allocate array of processor info blocks */

<<<<<<< HEAD
    pMPI_ProcessorInfo = tMPI_Malloc( sizeof(MPI_NUMA_PROCESSOR_INFO) * 
                                      dwTotalProcessors );
    if(pMPI_ProcessorInfo == NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"tMPI_Malloc failed for processor information");
=======
    pMPI_ProcessorInfo = tMPI_Malloc( sizeof(MPI_NUMA_PROCESSOR_INFO) *
                                      dwTotalProcessors );
    if (pMPI_ProcessorInfo == NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS, "tMPI_Malloc failed for processor information");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        goto cleanup;
    }

    /* zero fill to cover reserved must be-zero fields */
    memset(pMPI_ProcessorInfo, 0, sizeof(MPI_NUMA_PROCESSOR_INFO) * dwTotalProcessors);

<<<<<<< HEAD
    /* loop through each processor group, and for each group, capture the 
       processor numbers and NUMA node information. */

    for(GroupIndex = 0 ; GroupIndex < wActiveGroupCount ; GroupIndex++)
    {
        DWORD dwGroupProcessorCount;
        BYTE ProcessorIndex;

        dwGroupProcessorCount = func_GetActiveProcessorCount( GroupIndex );

        for(ProcessorIndex = 0 ; ProcessorIndex < dwGroupProcessorCount ; 
            ProcessorIndex++)
        {
            PROCESSOR_NUMBER *pProcessorNumber = &(pMPI_ProcessorInfo[i].ProcessorNumber);
            GROUP_AFFINITY *pGroupAffinity = &(pMPI_ProcessorInfo[i].GroupAffinity);
            USHORT *pNodeNumber = &(pMPI_ProcessorInfo[i].NumaNodeNumber);

            pProcessorNumber->Group = GroupIndex;
            pProcessorNumber->Number = ProcessorIndex;

            /* save an index to the processor array entry for the current processor
               this is used to enable subsequent threads to be created in a round 
               robin fashion starting at the next array entry
            */

            if( (CurrentProcessorNumber.Group == pProcessorNumber->Group ) &&
                (CurrentProcessorNumber.Number == pProcessorNumber->Number) )
=======
    /* loop through each processor group, and for each group, capture the
       processor numbers and NUMA node information. */

    for (GroupIndex = 0; GroupIndex < wActiveGroupCount; GroupIndex++)
    {
        DWORD dwGroupProcessorCount;
        BYTE  ProcessorIndex;

        dwGroupProcessorCount = func_GetActiveProcessorCount( GroupIndex );

        for (ProcessorIndex = 0; ProcessorIndex < dwGroupProcessorCount;
             ProcessorIndex++)
        {
            PROCESSOR_NUMBER *pProcessorNumber = &(pMPI_ProcessorInfo[i].ProcessorNumber);
            GROUP_AFFINITY   *pGroupAffinity   = &(pMPI_ProcessorInfo[i].GroupAffinity);
            USHORT           *pNodeNumber      = &(pMPI_ProcessorInfo[i].NumaNodeNumber);

            pProcessorNumber->Group  = GroupIndex;
            pProcessorNumber->Number = ProcessorIndex;

            /* save an index to the processor array entry for the current processor
               this is used to enable subsequent threads to be created in a round
               robin fashion starting at the next array entry
             */

            if ( (CurrentProcessorNumber.Group == pProcessorNumber->Group ) &&
                 (CurrentProcessorNumber.Number == pProcessorNumber->Number) )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                /* set global: current thread index into processor array */
                g_ulThreadIndex = i;
            }

            /* capture the node number and group affinity associated with processor entry
<<<<<<< HEAD
               any failures here are assumed to be catastrophic and disable 
               the group & NUMA aware thread support
            */

            if(!func_GetNumaProcessorNodeEx(pProcessorNumber, pNodeNumber))
=======
               any failures here are assumed to be catastrophic and disable
               the group & NUMA aware thread support
             */

            if (!func_GetNumaProcessorNodeEx(pProcessorNumber, pNodeNumber))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                tMPI_Fatal_error(TMPI_FARGS,
                                 "Processor enumeration, GetNumaProcessorNodeEx failed, error code=%d",
                                 GetLastError());
                goto cleanup;
            }

<<<<<<< HEAD
            if(!func_GetNumaNodeProcessorMaskEx(*pNodeNumber, pGroupAffinity))
=======
            if (!func_GetNumaNodeProcessorMaskEx(*pNodeNumber, pGroupAffinity))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                tMPI_Fatal_error(TMPI_FARGS,
                                 "Processor enumeration, GetNumaNodeProcessorMaskEx failed, error code=%d",
                                 GetLastError());
                goto cleanup;
            }

<<<<<<< HEAD
            /* future enhancement: construct GroupAffinity (single) processor 
=======
            /* future enhancement: construct GroupAffinity (single) processor
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
               mask within NUMA node for this processor entry */

            /* increment processor array index */
            i++;

            /* sanity check, should never happen */

<<<<<<< HEAD
            if(i > dwTotalProcessors)
            {
                tMPI_Fatal_error(TMPI_FARGS,"Processor enumeration exceeds allocated memory!");
=======
            if (i > dwTotalProcessors)
            {
                tMPI_Fatal_error(TMPI_FARGS, "Processor enumeration exceeds allocated memory!");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                goto cleanup;
            }
        }
    }

<<<<<<< HEAD
#if 0
    /* set the NUMA node affinity for the current thread
       failures to set the current thread affinity are ignored, 
       as a fringe case can arise on >32 processor systems with a 32bit 
       build/code.
    */
    func_SetThreadIdealProcessorEx(GetCurrentThread(), 
                                   &CurrentProcessorNumber, 
                                   NULL);

    if(func_GetNumaProcessorNodeEx(&CurrentProcessorNumber, 
                                   &CurrentNumaNodeNumber))
    {
        /* for the NUMA node number associated with the current processor 
           number, get the group affinity mask */
        if(func_GetNumaNodeProcessorMaskEx(CurrentNumaNodeNumber, 
                                           &CurrentThreadGroupAffinity))
        {
            /* set the current thread affinity to prevent it from running on 
               other NUMA nodes */
            func_SetThreadGroupAffinity(GetCurrentThread(), 
                                        &CurrentThreadGroupAffinity, 
                                        NULL);
        }
    }
#endif
 
    /* capture number of processors, highest NUMA node number, and processor 
       array */
    g_ulTotalProcessors = dwTotalProcessors;
    g_ulHighestNumaNodeNumber = ulHighestNumaNodeNumber;
    g_MPI_ProcessorInfo = pMPI_ProcessorInfo;

    iRet = 0 ;

#if 0   
    // TODO: debug DISCARD                        
    printf("primary thread tid=%lu group=%lu mask=0x%I64x group=%lu number=%lu ulThreadIndex=%lu\n",
        GetCurrentThreadId(),
        CurrentThreadGroupAffinity.Group,
        (ULONGLONG)CurrentThreadGroupAffinity.Mask,
        (ULONG)CurrentProcessorNumber.Group,
        (ULONG)CurrentProcessorNumber.Number,
        g_ulThreadIndex);
#endif

cleanup:

    if( iRet != 0 )
    {
        if( pMPI_ProcessorInfo )
=======

    /* capture number of processors, highest NUMA node number, and processor
       array */
    g_ulTotalProcessors       = dwTotalProcessors;
    g_ulHighestNumaNodeNumber = ulHighestNumaNodeNumber;
    g_MPI_ProcessorInfo       = pMPI_ProcessorInfo;

    iRet = 0;

cleanup:

    if (iRet != 0)
    {
        if (pMPI_ProcessorInfo)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            tMPI_Free( pMPI_ProcessorInfo );
        }
    }

    return 0;
}

<<<<<<< HEAD
=======
static void tMPI_Thread_id_list_init(void)
{
    EnterCriticalSection( &thread_id_list_lock );

    N_thread_id_list      = 0;
    Nalloc_thread_id_list = 4; /* number of initial allocation*/
    thread_id_list        = (thread_id_list_t*)tMPI_Malloc(
                sizeof(thread_id_list_t)*
                Nalloc_thread_id_list);

    LeaveCriticalSection( &thread_id_list_lock );
}


/* add an entry to the thread ID list, assuming it's locked */
static void tMPI_Thread_id_list_add_locked(DWORD               thread_id,
                                           struct tMPI_Thread *th)
{
    if (Nalloc_thread_id_list < N_thread_id_list + 1)
    {
        thread_id_list_t* new_list;
        int               i;

        /* double the size */
        Nalloc_thread_id_list *= 2;
        new_list               = (thread_id_list_t*)tMPI_Malloc(
                    sizeof(thread_id_list_t)*
                    Nalloc_thread_id_list);
        /* and copy over all elements */
        for (i = 0; i < N_thread_id_list; i++)
        {
            new_list[i] = thread_id_list[i];
        }
        /* free the old list */
        tMPI_Free(thread_id_list);
        thread_id_list = new_list;
    }
    thread_id_list[ N_thread_id_list ].thread_id = thread_id;
    thread_id_list[ N_thread_id_list ].th        = th;
    N_thread_id_list++;


}


/* add an entry to the thread ID list */
static void tMPI_Thread_id_list_add(DWORD thread_id, struct tMPI_Thread *th)
{
    EnterCriticalSection( &thread_id_list_lock );
    tMPI_Thread_id_list_add_locked(thread_id, th);
    LeaveCriticalSection( &thread_id_list_lock );
}

/* Remove an entry from the thread_id list, assuming it's locked */
static void tMPI_Thread_id_list_remove_locked(DWORD thread_id)
{
    int       i;
    tmpi_bool found = FALSE;

    /* move the last thread_id_list item to the one we want to remove */
    for (i = 0; i < N_thread_id_list; i++)
    {
        if (thread_id_list[i].thread_id == thread_id)
        {
            thread_id_list[i] = thread_id_list[N_thread_id_list - 1];
            found             = TRUE;
            break;
        }
    }

    if (found)
    {
        N_thread_id_list--;
    }
}


/* Remove an entry from the thread_id list */
static void tMPI_Thread_id_list_remove(DWORD thread_id)
{

    EnterCriticalSection( &thread_id_list_lock );
    tMPI_Thread_id_list_remove_locked(thread_id);
    LeaveCriticalSection( &thread_id_list_lock );
}



/* try to find a thread id in the thread id list. Return NULL when there is no
   such thread id in the list. Assumes the list is locked.*/
static struct tMPI_Thread *tMPI_Thread_id_list_find_locked(DWORD thread_id)
{
    int                 i;
    struct tMPI_Thread *ret = NULL;

    /* this is a linear search but it's only O(Nthreads). */
    for (i = 0; i < N_thread_id_list; i++)
    {
        if (thread_id_list[i].thread_id == thread_id)
        {
            ret = thread_id_list[i].th;
            break;
        }
    }

    return ret;
}

/* try to find a thread id in the thread id list. Return NULL when there is no
   such thread id in the list.*/
static struct tMPI_Thread *tMPI_Thread_id_list_find(DWORD thread_id)
{
    struct tMPI_Thread *ret = NULL;

    EnterCriticalSection( &thread_id_list_lock );
    ret = tMPI_Thread_id_list_find_locked(thread_id);

    LeaveCriticalSection( &thread_id_list_lock );
    return ret;
}

/* try to add the running thread to the list. Returns the tMPI_Thrread struct
   associated with this thread.*/
static struct tMPI_Thread *tMPI_Thread_id_list_add_self(void)
{
    DWORD               thread_id;
    struct tMPI_Thread *th = NULL;

    EnterCriticalSection( &thread_id_list_lock );

    thread_id = GetCurrentThreadId();
    th        = tMPI_Thread_id_list_find_locked(thread_id);
    if (th == NULL)
    {
        /* if not, create an ID, set it and return it */
        th = (struct tMPI_Thread*)tMPI_Malloc(sizeof(struct tMPI_Thread)*1);

        /* to create a handle that can be used outside of the current
           thread, the handle from GetCurrentThread() must first
           be duplicated.. */
        DuplicateHandle(GetCurrentProcess(),
                        GetCurrentThread(),
                        GetCurrentProcess(),
                        &th->th,
                        0,
                        FALSE,
                        DUPLICATE_SAME_ACCESS);

        /* This causes a small memory leak that is hard to fix. */
        th->started_by_tmpi = 0;
        tMPI_Thread_id_list_add_locked(thread_id, th);
    }
    LeaveCriticalSection( &thread_id_list_lock );

    return th;
}


>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
static void tMPI_Init_initers(void)
{
    int state;
    /* we can pre-check because it's atomic */
    if (tMPI_Atomic_get(&init_inited) == 0)
    {
        /* this can be a spinlock because the chances of collision are low. */
        tMPI_Spinlock_lock( &init_init );

<<<<<<< HEAD
        state=tMPI_Atomic_get(&init_inited);
=======
        state = tMPI_Atomic_get(&init_inited);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        tMPI_Atomic_memory_barrier_acq();
        if (state == 0)
        {
            InitializeCriticalSection(&mutex_init);
            InitializeCriticalSection(&once_init);
            InitializeCriticalSection(&cond_init);
            InitializeCriticalSection(&barrier_init);
<<<<<<< HEAD

            /* fatal errors are handled by the routine by calling tMPI_Fatal_error() */
            tMPI_Init_NUMA();	
=======
            InitializeCriticalSection(&thread_id_list_lock);

            /* fatal errors are handled by the routine by calling
               tMPI_Fatal_error() */
            tMPI_Init_NUMA();

            tMPI_Thread_id_list_init();
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

            tMPI_Atomic_memory_barrier_rel();
            tMPI_Atomic_set(&init_inited, 1);
        }

        tMPI_Spinlock_unlock( &init_init );
    }
}

<<<<<<< HEAD
=======


>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* TODO: this needs to go away!  (there's another one in pthreads.c)
   fatal errors are thankfully really rare*/
void tMPI_Fatal_error(const char *file, int line, const char *message, ...)
{
    va_list ap;

    fprintf(stderr, "tMPI Fatal error in %s, line %d: ", file, line);
    va_start(ap, message);
    vfprintf(stderr, message, ap);
    va_end(ap);
<<<<<<< HEAD
    fprintf(stderr,"\n");
=======
    fprintf(stderr, "\n");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    abort();
}



enum tMPI_Thread_support tMPI_Thread_support(void)
{
    return TMPI_THREAD_SUPPORT_YES;
}

struct tMPI_Thread_starter_param
{
<<<<<<< HEAD
    void *(*start_routine)(void*); /* the function */
    void *param; /* its parameter */
};

static DWORD WINAPI tMPI_Win32_thread_starter( LPVOID lpParam ) 
{
    struct tMPI_Thread_starter_param *prm=
              (struct tMPI_Thread_starter_param*)lpParam;
=======
    void               *(*start_routine)(void*); /* the function */
    void               *param;                   /* its parameter */
    struct tMPI_Thread *thread;
};

static DWORD WINAPI tMPI_Win32_thread_starter( LPVOID lpParam )
{
    struct tMPI_Thread_starter_param *prm =
        (struct tMPI_Thread_starter_param*)lpParam;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    (prm->start_routine)(prm->param);
    return 0;
}

<<<<<<< HEAD
HANDLE tMPI_Thread_create_NUMA(LPSECURITY_ATTRIBUTES lpThreadAttributes,
                               SIZE_T dwStackSize,
                               LPTHREAD_START_ROUTINE lpStartAddress,
                               LPVOID lpParameter,
                               DWORD dwCreationFlags,
                               LPDWORD lpThreadId)
{
    LPPROC_THREAD_ATTRIBUTE_LIST pAttributeList = NULL;
    HANDLE hThread = NULL;
    SIZE_T cbAttributeList = 0;
    GROUP_AFFINITY GroupAffinity;
    PROCESSOR_NUMBER IdealProcessorNumber;
    ULONG CurrentProcessorIndex;

    /* for each thread created, round-robin through the set of valid 
       processors and affinity masks.
       the assumption is that callers of tMPI_Thread_create_NUMA are creating 
       threads that saturate a given processor.
       for cases where threads are being created that rarely do work, standard 
       thread creation (eg: CreateThread) should be invoked instead.
    */

    CurrentProcessorIndex = (ULONG)InterlockedIncrement((volatile LONG *)&g_ulThreadIndex);
    CurrentProcessorIndex = CurrentProcessorIndex % g_ulTotalProcessors;

    /* group, mask. */

    memcpy(&GroupAffinity, 
           &(g_MPI_ProcessorInfo[CurrentProcessorIndex].GroupAffinity), 
           sizeof(GROUP_AFFINITY));

    /* group, processor number */
    
    memcpy(&IdealProcessorNumber, 
           &(g_MPI_ProcessorInfo[CurrentProcessorIndex].ProcessorNumber), 
           sizeof(PROCESSOR_NUMBER)); 

    /* determine size of allocation for AttributeList */

    if(!func_InitializeProcThreadAttributeList(pAttributeList,
                                               2,
                                               0,
                                               &cbAttributeList))
    {
        DWORD dwLastError = GetLastError();
        if( dwLastError != ERROR_INSUFFICIENT_BUFFER )
        {
            tMPI_Fatal_error(TMPI_FARGS,
                             "InitializeProcThreadAttributeList, error code=%d",
                             dwLastError);
            goto cleanup;
        }
    }

    pAttributeList = (LPPROC_THREAD_ATTRIBUTE_LIST)tMPI_Malloc( cbAttributeList );
    if( pAttributeList == NULL )
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to allocate pAttributeList");
        goto cleanup;
    }

    memset( pAttributeList, 0, cbAttributeList );

    if(!func_InitializeProcThreadAttributeList(pAttributeList,
                                               2,
                                               0,
                                               &cbAttributeList))
    {
        tMPI_Fatal_error(TMPI_FARGS,
                         "InitializeProcThreadAttributeList, error code=%d",
                         GetLastError());
        goto cleanup;
    }

    if(!func_UpdateProcThreadAttribute(pAttributeList,
                                       0,
                                       PROC_THREAD_ATTRIBUTE_GROUP_AFFINITY,
                                       &GroupAffinity,
                                       sizeof(GroupAffinity),
                                       NULL,
                                       NULL))
    {
        tMPI_Fatal_error(TMPI_FARGS,"UpdateProcThreadAttribute, error code=%d",
                         GetLastError());
        goto cleanup;
    }

    if(!func_UpdateProcThreadAttribute(pAttributeList,
                                       0,
                                       PROC_THREAD_ATTRIBUTE_IDEAL_PROCESSOR,
                                       &IdealProcessorNumber,
                                       sizeof(IdealProcessorNumber),
                                       NULL,
                                       NULL))
    {
        tMPI_Fatal_error(TMPI_FARGS,"UpdateProcThreadAttribute, error code=%d",
                         GetLastError());
        goto cleanup;
    }


    hThread = func_CreateRemoteThreadEx( GetCurrentProcess(),
                                         lpThreadAttributes,
                                         dwStackSize,
                                         lpStartAddress,
                                         lpParameter,
                                         dwCreationFlags,
                                         pAttributeList,
                                         lpThreadId);
            
    func_DeleteProcThreadAttributeList( pAttributeList );

#if 0   
	// TODO: debug only or DISCARD
    if( hThread )
    {
        PROCESSOR_NUMBER ProcNumber;
        USHORT NodeNumber;

        GetThreadIdealProcessorEx(hThread, &ProcNumber);
        GetNumaProcessorNodeEx(&ProcNumber, &NodeNumber);

        printf("started thread tid=%lu group=%lu mask=0x%I64x number=%lu numanode=%lu\n",
            *lpThreadId,
            GroupAffinity.Group,
            (ULONGLONG)GroupAffinity.Mask,
            ProcNumber.Number,
            NodeNumber
            );
    }
#endif

cleanup:
    
    if( pAttributeList )
    {
        tMPI_Free( pAttributeList );
    }

    return hThread;
}

int tMPI_Thread_get_hw_number(void)
{
    int ret;
=======

int tMPI_Thread_get_hw_number(void)
{
    int         ret;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    SYSTEM_INFO sysinfo;
    GetSystemInfo( &sysinfo );

<<<<<<< HEAD
    ret=sysinfo.dwNumberOfProcessors;
=======
    ret = sysinfo.dwNumberOfProcessors;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return ret;
}


<<<<<<< HEAD
=======


>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
int tMPI_Thread_create(tMPI_Thread_t *thread,
                       void *(*start_routine)(void *), void *arg)
{
    DWORD thread_id;
    struct tMPI_Thread_starter_param *prm;

    tMPI_Init_initers();

<<<<<<< HEAD
    /* a small memory leak to be sure that it doesn't get deallocated 
       once this function ends, before the newly created thread uses it. */
    prm=(struct tMPI_Thread_starter_param*)
              tMPI_Malloc(sizeof(struct tMPI_Thread_starter_param));
    prm->start_routine= start_routine;
    prm->param=arg;

    *thread=(struct tMPI_Thread*)tMPI_Malloc(sizeof(struct tMPI_Thread)*1);

    if(thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid thread pointer.");
        return EINVAL;
    }
    /* just create a plain thread. */
    (*thread)->th = CreateThread(NULL,
                                 0,
                                 tMPI_Win32_thread_starter,
                                 prm,
                                 0, 
                                 &thread_id);

    if((*thread)->th==NULL)
    {
        tMPI_Free(thread);
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create thread, error code=%d",
                         GetLastError());
        return -1;
    }

    /* inherit the thread priority from the parent thread. */
    /* TODO: is there value in setting this, vs. just allowing it to default 
       from the process?  currently, this limits the effectivenes of changing 
=======
    /* a small memory leak to be sure that it doesn't get deallocated
       once this function ends, before the newly created thread uses it. */
    prm = (struct tMPI_Thread_starter_param*)
        tMPI_Malloc(sizeof(struct tMPI_Thread_starter_param));
    prm->start_routine = start_routine;
    prm->param         = arg;

    *thread = (struct tMPI_Thread*)tMPI_Malloc(sizeof(struct tMPI_Thread)*1);

    if (thread == NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS, "Invalid thread pointer.");
        return EINVAL;
    }
    /* this must be locked before the thread is created to prevent a race
       condition if the thread immediately wants to create its own entry */
    EnterCriticalSection( &thread_id_list_lock );
    /* just create a plain thread. */
    (*thread)->started_by_tmpi = 1;
    (*thread)->th              = CreateThread(NULL,
                                              0,
                                              tMPI_Win32_thread_starter,
                                              prm,
                                              0,
                                              &thread_id);
    (*thread)->id = thread_id;

    if ((*thread)->th == NULL)
    {
        tMPI_Free(thread);
        tMPI_Fatal_error(TMPI_FARGS, "Failed to create thread, error code=%d",
                         GetLastError());
        return -1;
    }
    tMPI_Thread_id_list_add_locked(thread_id, (*thread));
    LeaveCriticalSection( &thread_id_list_lock );

    /* inherit the thread priority from the parent thread. */
    /* TODO: is there value in setting this, vs. just allowing it to default
       from the process?  currently, this limits the effectivenes of changing
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
       the priority in eg: TaskManager. */
    SetThreadPriority(((*thread)->th), GetThreadPriority(GetCurrentThread()));

    return 0;
}





<<<<<<< HEAD
int tMPI_Thread_create_aff(tMPI_Thread_t *thread,
                           void *(*start_routine)(void *), void *arg)
{
    DWORD thread_id;
    struct tMPI_Thread_starter_param *prm;

    tMPI_Init_initers();
    tMPI_Set_main_thread_affinity();

    /* a small memory leak to be sure that it doesn't get deallocated 
       once this function ends, before the newly created thread uses it. */
    prm=(struct tMPI_Thread_starter_param*)
              tMPI_Malloc(sizeof(struct tMPI_Thread_starter_param));
    prm->start_routine= start_routine;
    prm->param=arg;

    *thread=(struct tMPI_Thread*)tMPI_Malloc(sizeof(struct tMPI_Thread)*1);

    if(thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid thread pointer.");
        return EINVAL;
    }

    if( g_ulHighestNumaNodeNumber != 0 )
    {
        /* if running on a NUMA system, use the group and NUMA aware thread 
           creation logic */
        (*thread)->th = tMPI_Thread_create_NUMA(NULL,
                                                0,
                                                tMPI_Win32_thread_starter,
                                                prm,
                                                0, 
                                                &thread_id);
    } else {
        /* TODO: for now, non-NUMA systems don't set thread affinity. */
        (*thread)->th = CreateThread(NULL,
                                     0,
                                     tMPI_Win32_thread_starter,
                                     prm,
                                     0, 
                                     &thread_id);
    }

    if((*thread)->th==NULL)
    {
        tMPI_Free(thread);
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create thread, error code=%d",
                         GetLastError());
        return -1;
    }

    /* inherit the thread priority from the parent thread. */
    /* TODO: is there value in setting this, vs. just allowing it to default 
       from the process?  currently, this limits the effectivenes of changing 
       the priority in eg: TaskManager. */
    SetThreadPriority(((*thread)->th), GetThreadPriority(GetCurrentThread()));

    return 0;
}

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


int tMPI_Thread_join(tMPI_Thread_t thread, void **value_ptr)
{
<<<<<<< HEAD
    DWORD ret,retval;
=======
    DWORD ret, retval;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    ret = WaitForSingleObject(thread->th, INFINITE);

    if (ret != 0)
    {
<<<<<<< HEAD
        tMPI_Fatal_error(TMPI_FARGS,"Failed to join thread. error code=%d",
=======
        tMPI_Fatal_error(TMPI_FARGS, "Failed to join thread. error code=%d",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                         GetLastError());
        return -1;
    }

    if (value_ptr)
    {
        if (!GetExitCodeThread(thread, &retval))
        {
            /* TODO: somehow assign value_ptr */
            tMPI_Fatal_error(TMPI_FARGS,
                             "Failed to get thread exit code: error=%d",
                             GetLastError());
            return -1;
        }
    }
    CloseHandle(thread->th);
<<<<<<< HEAD
=======
    tMPI_Thread_id_list_remove(thread->id);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Free(thread);

    return 0;
}


void tMPI_Thread_exit(void *value_ptr)
{
    /* TODO: fix exit code */
    /* TODO: call destructors for thread-local storage */
    ExitThread( 0 );
}




int tMPI_Thread_cancel(tMPI_Thread_t thread)
{
    if (!TerminateThread( thread, -1) )
    {
<<<<<<< HEAD
        tMPI_Fatal_error(TMPI_FARGS,"Failed thread_cancel, error code=%d",
                         GetLastError());
        return -1;
    }
=======
        tMPI_Fatal_error(TMPI_FARGS, "Failed thread_cancel, error code=%d",
                         GetLastError());
        return -1;
    }
    tMPI_Thread_id_list_remove(thread->id);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return 0;
}


<<<<<<< HEAD


int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx) 
{
    if(mtx==NULL)
=======
tMPI_Thread_t tMPI_Thread_self(void)
{
    tMPI_Thread_t th;
    tMPI_Init_initers();

    th = tMPI_Thread_id_list_add_self();

    return th;
}


int tMPI_Thread_equal(tMPI_Thread_t t1, tMPI_Thread_t t2)
{
    /* because the tMPI thread IDs are unique, we can compare them directly */
    return (t1 == t2);
}

enum tMPI_Thread_setaffinity_support tMPI_Thread_setaffinity_support(void)
{
    /* Windows supports seting of thread affinities */
    return TMPI_SETAFFINITY_SUPPORT_YES;
}

int tMPI_Thread_setaffinity_single(tMPI_Thread_t thread, unsigned int nr)
{
    GROUP_AFFINITY   GroupAffinity;
    PROCESSOR_NUMBER IdealProcessorNumber;
    /* thread NUMA node */
    USHORT           NumaNodeNumber;

    /* check for a processor info array. This exists if NUMA
       style calls have been succesfully initialized. */
    if (g_MPI_ProcessorInfo != NULL)
    {

        /*func_GetCurrentProcessorNumberEx(&CurrentProcessorNumber);*/
        /* group, mask. */
        memcpy(&GroupAffinity,
               &(g_MPI_ProcessorInfo[nr].GroupAffinity),
               sizeof(GROUP_AFFINITY));

        /* group, processor number */

        memcpy(&IdealProcessorNumber,
               &(g_MPI_ProcessorInfo[nr].ProcessorNumber),
               sizeof(PROCESSOR_NUMBER));


        /* set the NUMA node affinity for the current thread
           failures to set the current thread affinity are ignored,
           as a fringe case can arise on >32 processor systems with a 32bit
           build/code.
         */
        func_SetThreadIdealProcessorEx(thread->th,
                                       &IdealProcessorNumber,
                                       NULL);

        if (func_GetNumaProcessorNodeEx(&IdealProcessorNumber,
                                        &NumaNodeNumber))
        {
            /* for the NUMA node number associated with the current processor
               number, get the group affinity mask */
            if (func_GetNumaNodeProcessorMaskEx(NumaNodeNumber,
                                                &GroupAffinity))
            {
                /* set the current thread affinity to prevent it from running
                   on other NUMA nodes */
                func_SetThreadGroupAffinity(thread->th,
                                            &GroupAffinity,
                                            NULL);
                return 0;
            }
        }
        return 1;
    }
    else
    {
        /* No NUMA-style calls. We just do a simpler thing. */
        if ( (func_SetThreadIdealProcessor != NULL) )
        {
            return (func_SetThreadIdealProcessor(thread->th, nr) == -1);
        }
    }
    return 0;
}



int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx)
{
    if (mtx == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return EINVAL;
    }

<<<<<<< HEAD
    mtx->mutex=(struct tMPI_Mutex*)tMPI_Malloc(sizeof(struct tMPI_Mutex)*1);
=======
    mtx->mutex = (struct tMPI_Mutex*)tMPI_Malloc(sizeof(struct tMPI_Mutex)*1);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    InitializeCriticalSection(&(mtx->mutex->cs));

    return 0;
}


<<<<<<< HEAD
int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx) 
{
    if(mtx == NULL)
=======
int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx)
{
    if (mtx == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return EINVAL;
    }

    DeleteCriticalSection(&(mtx->mutex->cs));
    tMPI_Free(mtx->mutex);

    return 0;
}




static int tMPI_Thread_mutex_init_once(tMPI_Thread_mutex_t *mtx)
{
<<<<<<< HEAD
    int ret=0;
=======
    int ret = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the mutex init routine instead.
     * It might seem like overkill, but it will only be executed the first
<<<<<<< HEAD
     * time you call a static mutex, and it is important to get all the 
     * memory barriers right. Trust me, you don't want a deadlock here...
     */ 
=======
     * time you call a static mutex, and it is important to get all the
     * memory barriers right. Trust me, you don't want a deadlock here...
     */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* initialize the initializers */
    tMPI_Init_initers();
    /* Lock the common one-time init mutex so we can check carefully */
    EnterCriticalSection( &mutex_init );

    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (mtx->mutex == NULL)
    {
        /* No need to keep the lock during execution -
         * Only one thread can do it anyway.
         */
<<<<<<< HEAD
        ret=tMPI_Thread_mutex_init(mtx);
=======
        ret = tMPI_Thread_mutex_init(mtx);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    LeaveCriticalSection( &mutex_init );

    return ret;
}



int tMPI_Thread_mutex_lock(tMPI_Thread_mutex_t *mtx)
{
    /* check whether the mutex is initialized */
    if (tMPI_Atomic_get( &(mtx->initialized)  ) == 0)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }

    /* The mutex is now guaranteed to be valid. */
    EnterCriticalSection( &(mtx->mutex->cs) );

    return 0;
}




int tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx)
{
    BOOL ret;

    /* check whether the mutex is initialized */
    if (tMPI_Atomic_get( &(mtx->initialized)  ) == 0)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }

    /* The mutex is now guaranteed to be valid. */
<<<<<<< HEAD
    ret=TryEnterCriticalSection( &(mtx->mutex->cs) );
=======
    ret = TryEnterCriticalSection( &(mtx->mutex->cs) );
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return (ret != 0);
}



int tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx)
{
    /* we should have initialized our critical section anyway */
    LeaveCriticalSection( &(mtx->mutex->cs) );

    return 0;
}



int tMPI_Thread_key_create(tMPI_Thread_key_t *key, void (*destructor)(void *))
{
<<<<<<< HEAD
    if(key==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid key pointer.");
=======
    if (key == NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS, "Invalid key pointer.");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        return EINVAL;
    }


    /* TODO: make list of destructors for thread-local storage */
<<<<<<< HEAD
    key->key=(struct tMPI_Thread_key*)tMPI_Malloc(sizeof(struct 
                                                         tMPI_Thread_key)*1);
 
    (key)->key->wkey=TlsAlloc();

    if ( (key)->key->wkey == TLS_OUT_OF_INDEXES ) 
=======
    key->key = (struct tMPI_Thread_key*)tMPI_Malloc(sizeof(struct
                                                           tMPI_Thread_key)*1);

    (key)->key->wkey = TlsAlloc();

    if ( (key)->key->wkey == TLS_OUT_OF_INDEXES)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        tMPI_Fatal_error(TMPI_FARGS,
                         "Failed to create thread key, error code=%d.",
                         GetLastError());
        return -1;
    }

    return 0;
}


int tMPI_Thread_key_delete(tMPI_Thread_key_t key)
{
    TlsFree(key.key->wkey);
    tMPI_Free(key.key);

    return 0;
}



void * tMPI_Thread_getspecific(tMPI_Thread_key_t key)
{
    void *p = NULL;

<<<<<<< HEAD
    p=TlsGetValue(key.key->wkey);
=======
    p = TlsGetValue(key.key->wkey);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return p;
}


int tMPI_Thread_setspecific(tMPI_Thread_key_t key, void *value)
{
    BOOL ret;

    ret = TlsSetValue(key.key->wkey, value);

<<<<<<< HEAD
    return ret==0;
=======
    return ret == 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

#if 0
/* use once Vista is minimum required version */
static BOOL CALLBACK InitHandleWrapperFunction(PINIT_ONCE InitOnce,
<<<<<<< HEAD
                                               PVOID Parameter,
                                               PVOID *lpContext)
{
    void (*fn)(void)=(void (*)(void))Parameter;
=======
                                               PVOID      Parameter,
                                               PVOID     *lpContext)
{
    void (*fn)(void) = (void (*)(void))Parameter;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fn();

    return TRUE;
}

CRITICAL_SECTION tMPI_Once_cs;
<<<<<<< HEAD
tMPI_Spinlock_t tMPI_Once_cs_lock=TMPI_SPINLOCK_INITIALIZER;
volatile int tMPI_Once_init=0;
#endif

int tMPI_Thread_once(tMPI_Thread_once_t *once_control, 
                     void (*init_routine)(void))
=======
tMPI_Spinlock_t  tMPI_Once_cs_lock = TMPI_SPINLOCK_INITIALIZER;
volatile int     tMPI_Once_init    = 0;
#endif

int tMPI_Thread_once(tMPI_Thread_once_t *once_control,
                     void                (*init_routine)(void))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
#if 0
    /* use once Vista is minimum required version */
    BOOL bStatus;
<<<<<<< HEAD
    bStatus = InitOnceExecuteOnce(once_control, InitHandleWrapperFunction, 
=======
    bStatus = InitOnceExecuteOnce(once_control, InitHandleWrapperFunction,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                  init_routine, NULL);

    if (!bStatus)
    {
<<<<<<< HEAD
        tMPI_Fatal_error(TMPI_FARGS,"Failed to run thread_once routine");
=======
        tMPI_Fatal_error(TMPI_FARGS, "Failed to run thread_once routine");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        return -1;
    }
#else
    /* really ugly hack - and it's slow... */
    tMPI_Init_initers();
    EnterCriticalSection(&once_init);
    if (tMPI_Atomic_get(&(once_control->once)) == 0)
    {
        (*init_routine)();
        tMPI_Atomic_set(&(once_control->once), 1);
    }
    LeaveCriticalSection(&once_init);
#endif
    return 0;
}





<<<<<<< HEAD
int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond) 
{
    if(cond==NULL)
=======
int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond)
{
    if (cond == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return EINVAL;
    }

<<<<<<< HEAD
    cond->condp=(struct tMPI_Thread_cond*)
              tMPI_Malloc(sizeof(struct tMPI_Thread_cond)*1);
=======
    cond->condp = (struct tMPI_Thread_cond*)
        tMPI_Malloc(sizeof(struct tMPI_Thread_cond)*1);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#if 0
    /* use this code once Vista is the minimum version required */
    InitializeConditionVariable( &(cond->cv) );
#else
<<<<<<< HEAD
    cond->condp->Nwaiters=0;
    InitializeCriticalSection(&(cond->condp->wtr_lock));
    cond->condp->Nrelease=0;
    cond->condp->cycle=0;
    /* a manual reset, unsignalled event */
    cond->condp->ev = CreateEvent(NULL, TRUE, FALSE, NULL); 
=======
    cond->condp->Nwaiters = 0;
    InitializeCriticalSection(&(cond->condp->wtr_lock));
    cond->condp->Nrelease = 0;
    cond->condp->cycle    = 0;
    /* a manual reset, unsignalled event */
    cond->condp->ev = CreateEvent(NULL, TRUE, FALSE, NULL);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    return 0;
}


<<<<<<< HEAD
int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond) 
=======
int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
#if 0
    /* use this code once Vista is the minimum version required */
    /* windows doesnt have this function */
#else
    DeleteCriticalSection(&(cond->condp->wtr_lock));
    tMPI_Free(cond->condp);
#endif
    return 0;
}



<<<<<<< HEAD
/*! \brief Static init routine for pthread barrier 
=======
/*! \brief Static init routine for pthread barrier
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for tMPI_Thread.h
<<<<<<< HEAD
 * 
 * \param cond  Condition variable, must be statically initialized
 *  
=======
 *
 * \param cond  Condition variable, must be statically initialized
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * \return status - 0 on success, or a standard error code.
 */
static int tMPI_Thread_cond_init_once(tMPI_Thread_cond_t *cond)
{
<<<<<<< HEAD
    int ret=0;
=======
    int ret = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the cond init routine instead.
     * It might seem like overkill, but it will only be executed the first
<<<<<<< HEAD
     * time you call a static condition variable, and it is important to get 
     * the memory barriers right. Trust me, you don't want a deadlock here...
     */ 
=======
     * time you call a static condition variable, and it is important to get
     * the memory barriers right. Trust me, you don't want a deadlock here...
     */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* initialize the initializers */
    tMPI_Init_initers();
    /* Lock the common one-time init mutex so we can check carefully */
    EnterCriticalSection( &cond_init );

    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (cond->condp == NULL)
    {
        /* No need to keep the lock during execution -
         * Only one thread can do it anyway.  */
<<<<<<< HEAD
        ret=tMPI_Thread_cond_init(cond);
=======
        ret = tMPI_Thread_cond_init(cond);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    LeaveCriticalSection( &cond_init );

    return ret;
}




int tMPI_Thread_cond_wait(tMPI_Thread_cond_t *cond, tMPI_Thread_mutex_t *mtx)
{
<<<<<<< HEAD
    BOOL wait_done=FALSE;
    BOOL last_waiter=FALSE;
    int my_cycle;
=======
    BOOL wait_done   = FALSE;
    BOOL last_waiter = FALSE;
    int  my_cycle;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        tMPI_Thread_cond_init_once(cond);
    }
    /* the mutex must have been initialized because it should be locked here */

#if 0
    /* use this code once Vista is the minimum version required */
<<<<<<< HEAD
    ret=SleepConditionVariableCS (&(cond->cv), &(mtx->cs), INFINITE);

    if (!ret)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed wait for condition, error code=%d",
=======
    ret = SleepConditionVariableCS (&(cond->cv), &(mtx->cs), INFINITE);

    if (!ret)
    {
        tMPI_Fatal_error(TMPI_FARGS, "Failed wait for condition, error code=%d",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                         GetLastError());
        return -1;
    }
#else
    /* serially increase waiter count */
    EnterCriticalSection(&(cond->condp->wtr_lock));
    cond->condp->Nwaiters++;
    my_cycle = cond->condp->cycle;
    LeaveCriticalSection(&(cond->condp->wtr_lock));

    /* now it's safe to release the mutex from the fn call */
    LeaveCriticalSection(&(mtx->mutex->cs));

    /* Loop a wait until we found out we've waited for the right event.
       Note that this loop is potentially a busy-wait loop in bad
       circumstances (higher priority threads, for example). */
    do
    {
        /* do the actual waiting */
<<<<<<< HEAD
        if (WaitForSingleObject( cond->condp->ev, INFINITE )== WAIT_FAILED)
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed event reset, error code=%d",
=======
        if (WaitForSingleObject( cond->condp->ev, INFINITE ) == WAIT_FAILED)
        {
            tMPI_Fatal_error(TMPI_FARGS, "Failed event reset, error code=%d",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                             GetLastError());
            return -1;
        }

        /* serially check whether we got the right event.  */
        EnterCriticalSection(&(cond->condp->wtr_lock));
<<<<<<< HEAD
        wait_done = (cond->condp->Nrelease > 0) && 
                    (cond->condp->cycle!=my_cycle);
        LeaveCriticalSection(&(cond->condp->wtr_lock));
    }
    while(!wait_done);
=======
        wait_done = (cond->condp->Nrelease > 0) &&
            (cond->condp->cycle != my_cycle);
        LeaveCriticalSection(&(cond->condp->wtr_lock));
    }
    while (!wait_done);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* We obtain the mutex from the function call */
    EnterCriticalSection(&(mtx->mutex->cs));

    /* we serially decrease the waiter count and release count */
    EnterCriticalSection(&(cond->condp->wtr_lock));
    cond->condp->Nwaiters--;
    cond->condp->Nrelease--;
<<<<<<< HEAD
    last_waiter=(cond->condp->Nrelease==0);
=======
    last_waiter = (cond->condp->Nrelease == 0);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    LeaveCriticalSection(&(cond->condp->wtr_lock));

    /* manually release the event if everybody's done with it */
    if (last_waiter)
    {
        if (!ResetEvent( cond->condp->ev ))
        {
<<<<<<< HEAD
            tMPI_Fatal_error(TMPI_FARGS,"Failed event reset, error code=%d",
=======
            tMPI_Fatal_error(TMPI_FARGS, "Failed event reset, error code=%d",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                             GetLastError());
            return -1;
        }
    }
#endif

    return 0;
}




int tMPI_Thread_cond_signal(tMPI_Thread_cond_t *cond)
{
    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        tMPI_Thread_cond_init_once(cond);
    }
    /* The condition variable is now guaranteed to be valid. */
#if 0
    /* use this code once Vista is the minimum version required */
    WakeConditionVariable( &(cond->cv) );
#else
    EnterCriticalSection(&(cond->condp->wtr_lock));
    /* check if we're not still busy with a release. If we are, do nothing. */
    if (cond->condp->Nwaiters > cond->condp->Nrelease)
    {
        cond->condp->Nrelease++;
        cond->condp->cycle++;
<<<<<<< HEAD
        if (!SetEvent(cond->condp->ev)) /* actually release the 
                                           waiting threads */
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed SetEvent, error code=%d",
=======
        if (!SetEvent(cond->condp->ev)) /* actually release the
                                           waiting threads */
        {
            tMPI_Fatal_error(TMPI_FARGS, "Failed SetEvent, error code=%d",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                             GetLastError());
            return -1;
        }
    }
    LeaveCriticalSection(&(cond->condp->wtr_lock));
#endif

    return 0;
}



int tMPI_Thread_cond_broadcast(tMPI_Thread_cond_t *cond)
{
    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        tMPI_Thread_cond_init_once(cond);
    }
    /* The condition variable is now guaranteed to be valid. */
#if 0
    /* use this code once Vista is the minimum version required */
    WakeAllConditionVariable( &(cond->cv) );
#else
    EnterCriticalSection(&(cond->condp->wtr_lock));
    /* check whether there are any waiters */
    if (cond->condp->Nwaiters > 0)
    {
<<<<<<< HEAD
        cond->condp->Nrelease=cond->condp->Nwaiters;
        cond->condp->cycle++;
        if (!SetEvent(cond->condp->ev)) /* actually release the 
                                           waiting threads */
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed SetEvent, error code=%d",
=======
        cond->condp->Nrelease = cond->condp->Nwaiters;
        cond->condp->cycle++;
        if (!SetEvent(cond->condp->ev)) /* actually release the
                                           waiting threads */
        {
            tMPI_Fatal_error(TMPI_FARGS, "Failed SetEvent, error code=%d",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                             GetLastError());
            return -1;
        }
    }
    LeaveCriticalSection(&(cond->condp->wtr_lock));
#endif
    return 0;
}




int tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *barrier, int n)
{
<<<<<<< HEAD
    if(barrier==NULL)
=======
    if (barrier == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return EINVAL;
    }

<<<<<<< HEAD
    barrier->barrierp=(struct tMPI_Thread_barrier*)
              tMPI_Malloc(sizeof(struct tMPI_Thread_barrier)*1);

#if 0
 /* use this once Vista is the oldest supported windows version: */
=======
    barrier->barrierp = (struct tMPI_Thread_barrier*)
        tMPI_Malloc(sizeof(struct tMPI_Thread_barrier)*1);

#if 0
    /* use this once Vista is the oldest supported windows version: */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    InitializeCriticalSection(&(barrier->barrierp->cs));
    InitializeConditionVariable(&(barrier->barrierp->cv));
#else
    tMPI_Thread_mutex_init(&(barrier->barrierp->cs));
    tMPI_Thread_cond_init(&(barrier->barrierp->cv));
#endif

    barrier->threshold = n;
    barrier->count     = n;
    barrier->cycle     = 0;

    return 0;
}



int tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *barrier)
<<<<<<< HEAD
{   
    if(barrier==NULL)
=======
{
    if (barrier == NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        return EINVAL;
    }

#if 0
    DeleteCriticalSection(&(barrier->barrierp->cs));
#else
    tMPI_Thread_mutex_destroy(&(barrier->barrierp->cs));
#endif

    tMPI_Thread_cond_destroy(&(barrier->barrierp->cv));

    tMPI_Free(barrier->barrierp);

    return 0;
}



<<<<<<< HEAD
/*! \brief Static init routine for pthread barrier 
=======
/*! \brief Static init routine for pthread barrier
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for tMPI_Thread.h
 *
 * \param barrier Statically initialized barrier type
 * \param n       Number of members in barrier
<<<<<<< HEAD
 * 
=======
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * \return status - 0 on success, or a standard error code.
 */
static int tMPI_Thread_barrier_init_once(tMPI_Thread_barrier_t *barrier, int n)
{
    int ret;

    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the cond init routine instead.
     * It might seem like overkill, but it will only be executed the first
<<<<<<< HEAD
     * time you call a static condition variable, and it is important to get 
     * the memory barriers right. Trust me, you don't want a deadlock here...
     */ 
=======
     * time you call a static condition variable, and it is important to get
     * the memory barriers right. Trust me, you don't want a deadlock here...
     */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    /* initialize the initializers */
    tMPI_Init_initers();

    /* Lock the common one-time init mutex so we can check carefully */
    EnterCriticalSection( &barrier_init );

    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (barrier->barrierp == NULL)
    {
        /* No need to keep the lock during execution -
         * Only one thread can do it anyway.  */
<<<<<<< HEAD
        ret=tMPI_Thread_barrier_init(barrier, n);
=======
        ret = tMPI_Thread_barrier_init(barrier, n);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    LeaveCriticalSection( &barrier_init );

    return ret;
}



int tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t *barrier)
{
<<<<<<< HEAD
    int    cycle;
    BOOL    rc=FALSE;
    int     ret=0;
=======
    int     cycle;
    BOOL    rc  = FALSE;
    int     ret = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /*tMPI_Thread_pthread_barrier_t *p;*/

    /* check whether the barrier is initialized */
    if (tMPI_Atomic_get( &(barrier->initialized)  ) == 0)
    {
<<<<<<< HEAD
        tMPI_Thread_barrier_init_once(barrier,barrier->threshold);        
=======
        tMPI_Thread_barrier_init_once(barrier, barrier->threshold);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

#if 0
    EnterCriticalSection( &(barrier->barrierp->cs)  );
#else
    tMPI_Thread_mutex_lock( &(barrier->barrierp->cs) );
#endif



    cycle = barrier->cycle;

    /* Decrement the count atomically and check if it is zero.
     * This will only be true for the last thread calling us.
     */
<<<<<<< HEAD
    if( --(barrier->count) <= 0 )
    { 
=======
    if (--(barrier->count) <= 0)
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        barrier->cycle = !barrier->cycle;
        barrier->count = barrier->threshold;
#if 0
        WakeAllConditionVariable( &(barrier->barrierp->cv) );
#else
        tMPI_Thread_cond_broadcast( &(barrier->barrierp->cv) );
#endif
    }
    else
    {
<<<<<<< HEAD
        while(cycle == barrier->cycle)
        {
#if 0
            rc=SleepConditionVariableCS (&(barrier->barrierp->cv), 
                                         &(barrier->barrierp->cs), 
                                         INFINITE);
            if(!rc) 
            {
                ret=-1;
=======
        while (cycle == barrier->cycle)
        {
#if 0
            rc = SleepConditionVariableCS (&(barrier->barrierp->cv),
                                           &(barrier->barrierp->cs),
                                           INFINITE);
            if (!rc)
            {
                ret = -1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                break;
            }
#else
            rc = tMPI_Thread_cond_wait(&barrier->barrierp->cv,
                                       &barrier->barrierp->cs);
<<<<<<< HEAD
            if(rc != 0) break;
=======
            if (rc != 0)
            {
                break;
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
        }
    }
#if 0
    LeaveCriticalSection( &(barrier->barrierp->cs)  );
#else
    tMPI_Thread_mutex_unlock( &(barrier->barrierp->cs) );
#endif
    return ret;
}

#else

/* just to have some symbols */
<<<<<<< HEAD
int tMPI_Thread_winthreads=0;

#endif /* THREAD_WINDOWS  */

=======
int tMPI_Thread_winthreads = 0;

#endif /* THREAD_WINDOWS  */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2