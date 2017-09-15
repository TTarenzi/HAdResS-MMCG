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

#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#ifdef HAVE_IO_H
#include <io.h>
#endif

#include "gmx_fatal.h"
#include "macros.h"
#include "smalloc.h"
#include "futil.h"
#include "filenm.h"
#include "string2.h"
#include "gmxfio.h"
#include "md5.h"

<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "thread_mpi.h"
#endif

#include "gmxfio_int.h"

/* This is the new improved and thread safe version of gmxfio. */



<<<<<<< HEAD
/* the list of open files is a linked list, with a dummy element at its head; 
=======
/* the list of open files is a linked list, with a dummy element at its head;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
       it is initialized when the first file is opened. */
static t_fileio *open_files = NULL;


<<<<<<< HEAD
#ifdef GMX_THREADS
/* this mutex locks the open_files structure so that no two threads can 
   modify it.        

   For now, we use this as a coarse grained lock on all file 
   insertion/deletion operations because it makes avoiding deadlocks 
=======
#ifdef GMX_THREAD_MPI
/* this mutex locks the open_files structure so that no two threads can
   modify it.

   For now, we use this as a coarse grained lock on all file
   insertion/deletion operations because it makes avoiding deadlocks
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   easier, and adds almost no overhead: the only overhead is during
   opening and closing of files, or during global operations like
   iterating along all open files. All these cases should be rare
   during the simulation. */
<<<<<<< HEAD
static tMPI_Thread_mutex_t open_file_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
=======
static tMPI_Thread_mutex_t open_file_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif


/* These simple lists define the I/O type for these files */
static const int ftpXDR[] =
<<<<<<< HEAD
    { efTPR, efTRR, efEDR, efXTC, efMTX, efCPT };
static const int ftpASC[] =
    { efTPA, efGRO, efPDB };
static const int ftpBIN[] =
    { efTPB, efTRJ };
#ifdef HAVE_XML
static const int ftpXML[] =
    {   efXML};
#endif

const char *itemstr[eitemNR] =
    { "[header]", "[inputrec]", "[box]", "[topology]", "[coordinates]",
                "[velocities]", "[forces]" };

const char *eioNames[eioNR] =
    { "REAL", "INT", "GMX_STE_T", "UCHAR", "NUCHAR", "USHORT", "RVEC", "NRVEC",
                "IVEC", "STRING" };
=======
{ efTPR, efTRR, efEDR, efXTC, efMTX, efCPT };
static const int ftpASC[] =
{ efTPA, efGRO, efPDB };
static const int ftpBIN[] =
{ efTPB, efTRJ };

const char *itemstr[eitemNR] =
{
    "[header]", "[inputrec]", "[box]", "[topology]", "[coordinates]",
    "[velocities]", "[forces]"
};

const char *eioNames[eioNR] =
{
    "REAL", "INT", "GMX_STE_T", "UCHAR", "NUCHAR", "USHORT", "RVEC", "NRVEC",
    "IVEC", "STRING"
};
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2



/* Comment strings for TPA only */
const char *comment_str[eitemNR] = {
<<<<<<< HEAD
"; The header holds information on the number of atoms etc. and on whether\n"
"; certain items are present in the file or not.\n"
"; \n"
";                             WARNING\n"
";                   DO NOT EDIT THIS FILE BY HAND\n"
"; The GROMACS preprocessor performs a lot of checks on your input that\n"
"; you ignore when editing this. Your simulation may crash because of this\n",
"; The inputrec holds the parameters for MD such as the number of steps,\n"
"; the timestep and the cut-offs.\n",
"; The simulation box in nm.\n",
"; The topology section describes the topology of the molecules\n"
"; i.e. bonds, angles and dihedrals etc. and also holds the force field\n"
"; parameters.\n", 
"; The atomic coordinates in nm\n",
"; The atomic velocities in nm/ps\n",
"; The forces on the atoms in nm/ps^2\n" };
=======
    "; The header holds information on the number of atoms etc. and on whether\n"
    "; certain items are present in the file or not.\n"
    "; \n"
    ";                             WARNING\n"
    ";                   DO NOT EDIT THIS FILE BY HAND\n"
    "; The GROMACS preprocessor performs a lot of checks on your input that\n"
    "; you ignore when editing this. Your simulation may crash because of this\n",
    "; The inputrec holds the parameters for MD such as the number of steps,\n"
    "; the timestep and the cut-offs.\n",
    "; The simulation box in nm.\n",
    "; The topology section describes the topology of the molecules\n"
    "; i.e. bonds, angles and dihedrals etc. and also holds the force field\n"
    "; parameters.\n",
    "; The atomic coordinates in nm\n",
    "; The atomic velocities in nm/ps\n",
    "; The forces on the atoms in nm/ps^2\n"
};
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2




/******************************************************************
 *
<<<<<<< HEAD
 * Internal functions: 
=======
 * Internal functions:
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *
 ******************************************************************/

static int gmx_fio_int_flush(t_fileio* fio)
{
    int rc = 0;

    if (fio->fp)
    {
        rc = fflush(fio->fp);
    }
    else if (fio->xdr)
    {
        rc = fflush((FILE *) fio->xdr->x_private);
    }

    return rc;
}

/* returns TRUE if the file type ftp is in the set set */
static gmx_bool in_ftpset(int ftp, int nset, const int set[])
{
<<<<<<< HEAD
    int i;
=======
    int      i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_bool bResult;

    bResult = FALSE;
    for (i = 0; (i < nset); i++)
<<<<<<< HEAD
        if (ftp == set[i])
            bResult = TRUE;
=======
    {
        if (ftp == set[i])
        {
            bResult = TRUE;
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return bResult;
}



extern void gmx_fio_set_comment(t_fileio *fio, const char *comment)
{
<<<<<<< HEAD
    fio->comment=comment;
=======
    fio->comment = comment;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

extern void gmx_fio_unset_comment(t_fileio *fio)
{
<<<<<<< HEAD
    fio->comment=NULL;
=======
    fio->comment = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


const char *gmx_fio_dbgstr(t_fileio *fio, const char *desc, char *buf)
{
    if (!fio->bDebug)
    {
        /* set to empty string */
<<<<<<< HEAD
        buf[0]=0;
=======
        buf[0] = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        snprintf(buf, GMX_FIO_BUFLEN, "  ; %s %s", fio->comment ? fio->comment : "", desc);
    }
    return buf;
}


/* check the number of items given against the type */
<<<<<<< HEAD
void gmx_fio_check_nitem(t_fileio *fio, int eio, int nitem, const char *file, 
                         int line)
{
    if ((nitem != 1) && !((eio == eioNRVEC) || (eio == eioNUCHAR)))
        gmx_fatal(FARGS, 
                  "nitem (%d) may differ from 1 only for %s or %s, not   for %s"
                  "(%s, %d)",nitem,eioNames[eioNUCHAR],eioNames[eioNRVEC],
                  eioNames[eio],file,line);
=======
void gmx_fio_check_nitem(t_fileio *fio, int eio, int nitem, const char *file,
                         int line)
{
    if ((nitem != 1) && !((eio == eioNRVEC) || (eio == eioNUCHAR)))
    {
        gmx_fatal(FARGS,
                  "nitem (%d) may differ from 1 only for %s or %s, not   for %s"
                  "(%s, %d)", nitem, eioNames[eioNUCHAR], eioNames[eioNRVEC],
                  eioNames[eio], file, line);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* output a data type error. */
<<<<<<< HEAD
void gmx_fio_fe(t_fileio *fio, int eio, const char *desc, 
=======
void gmx_fio_fe(t_fileio *fio, int eio, const char *desc,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                const char *srcfile, int line)
{

    gmx_fatal(FARGS, "Trying to %s %s type %d (%s), src %s, line %d",
<<<<<<< HEAD
              fio->bRead ? "read" : "write",desc,eio,
              ((eio >= 0) && (eio < eioNR)) ? eioNames[eio] : "unknown",
              srcfile,line);
=======
              fio->bRead ? "read" : "write", desc, eio,
              ((eio >= 0) && (eio < eioNR)) ? eioNames[eio] : "unknown",
              srcfile, line);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* set the reader/writer functions based on the file type */
static void gmx_fio_set_iotype(t_fileio *fio)
{
    if (in_ftpset(fio->iFTP, asize(ftpXDR), ftpXDR))
    {
<<<<<<< HEAD
#ifdef USE_XDR    
        fio->iotp=&xdr_iotype;
#else
        gmx_fatal(FARGS,"Sorry, no XDR");
=======
#ifdef USE_XDR
        fio->iotp = &xdr_iotype;
#else
        gmx_fatal(FARGS, "Sorry, no XDR");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    }
    else if (in_ftpset(fio->iFTP, asize(ftpASC), ftpASC))
    {
<<<<<<< HEAD
        fio->iotp=&asc_iotype;
    }
    else if (in_ftpset(fio->iFTP, asize(ftpBIN), ftpBIN))
    {
        fio->iotp=&bin_iotype;
    }
#ifdef HAVE_XMl
    else if (in_ftpset(fio->iFTP,asize(ftpXML),ftpXML))
    {
        fio->iotp=&dummy_iotype;
    }
#endif
    else
        fio->iotp=&dummy_iotype;
=======
        fio->iotp = &asc_iotype;
    }
    else if (in_ftpset(fio->iFTP, asize(ftpBIN), ftpBIN))
    {
        fio->iotp = &bin_iotype;
    }
    else
    {
        fio->iotp = &dummy_iotype;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}


/* lock the mutex associated with this fio. This needs to be done for every
   type of access to the fio's elements. */
void gmx_fio_lock(t_fileio *fio)
{
<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Lock_lock(&(fio->mtx));
#endif
}
/* unlock the mutex associated with this fio.  */
void gmx_fio_unlock(t_fileio *fio)
{
<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Lock_unlock(&(fio->mtx));
#endif
}

/* make a dummy head element, assuming we locked everything. */
static void gmx_fio_make_dummy(void)
{
    if (!open_files)
    {
<<<<<<< HEAD
        snew(open_files,1);
        open_files->fp=NULL;
        open_files->fn=NULL;
        open_files->next=open_files;
        open_files->prev=open_files;
#ifdef GMX_THREADS
=======
        snew(open_files, 1);
        open_files->fp   = NULL;
        open_files->fn   = NULL;
        open_files->next = open_files;
        open_files->prev = open_files;
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        tMPI_Lock_init(&(open_files->mtx));
#endif
    }
}







/***********************************************************************
 *
 * FILE LIST OPERATIONS
 *
<<<<<<< HEAD
***********************************************************************/
=======
 ***********************************************************************/
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/* insert a new t_fileio into the list */
static void gmx_fio_insert(t_fileio *fio)
{
    t_fileio *prev;
<<<<<<< HEAD
#ifdef GMX_THREADS
    /* first lock the big open_files mutex. */
    tMPI_Thread_mutex_lock(&open_file_mutex);
#endif
    /* now check whether the dummy element has been allocated, 
=======
#ifdef GMX_THREAD_MPI
    /* first lock the big open_files mutex. */
    tMPI_Thread_mutex_lock(&open_file_mutex);
#endif
    /* now check whether the dummy element has been allocated,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
       and allocate it if it hasn't */
    gmx_fio_make_dummy();

    /* and lock the fio we got and the list's head **/
    gmx_fio_lock(fio);
    gmx_fio_lock(open_files);
<<<<<<< HEAD
    prev=open_files->prev;
=======
    prev = open_files->prev;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* lock the element after the current one */
    if (prev != open_files)
    {
        gmx_fio_lock(prev);
    }

    /* now do the actual insertion: */
<<<<<<< HEAD
    fio->next=open_files;
    open_files->prev=fio;
    prev->next=fio;
    fio->prev=prev;
=======
    fio->next        = open_files;
    open_files->prev = fio;
    prev->next       = fio;
    fio->prev        = prev;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* now unlock all our locks */
    if (prev != open_files)
    {
        gmx_fio_unlock(prev);
    }
    gmx_fio_unlock(open_files);
    gmx_fio_unlock(fio);

<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* now unlock the big open_files mutex.  */
    tMPI_Thread_mutex_unlock(&open_file_mutex);
#endif
}

<<<<<<< HEAD
/* remove a t_fileio into the list. We assume the fio is locked, and we leave 
   it locked. 
   NOTE: We also assume that the open_file_mutex has been locked */
static void gmx_fio_remove(t_fileio *fio)
{    
    t_fileio *prev;

    /* lock prev, because we're changing it */ 
    gmx_fio_lock(fio->prev);

    /* now set the prev's pointer */
    fio->prev->next=fio->next;
=======
/* remove a t_fileio into the list. We assume the fio is locked, and we leave
   it locked.
   NOTE: We also assume that the open_file_mutex has been locked */
static void gmx_fio_remove(t_fileio *fio)
{
    t_fileio *prev;

    /* lock prev, because we're changing it */
    gmx_fio_lock(fio->prev);

    /* now set the prev's pointer */
    fio->prev->next = fio->next;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio->prev);

    /* with the next ptr, we can simply lock while the original was locked */
    gmx_fio_lock(fio->next);
<<<<<<< HEAD
    fio->next->prev=fio->prev;
    gmx_fio_unlock(fio->next);

    /* and make sure we point nowhere in particular */
    fio->next=fio->prev=fio;
}


/* get the first open file, or NULL if there is none. 
=======
    fio->next->prev = fio->prev;
    gmx_fio_unlock(fio->next);

    /* and make sure we point nowhere in particular */
    fio->next = fio->prev = fio;
}


/* get the first open file, or NULL if there is none.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   Returns a locked fio. */
static t_fileio *gmx_fio_get_first(void)
{
    t_fileio *ret;
    /* first lock the big open_files mutex and the dummy's mutex */

<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* first lock the big open_files mutex. */
    tMPI_Thread_mutex_lock(&open_file_mutex);
#endif
    gmx_fio_make_dummy();

    gmx_fio_lock(open_files);
<<<<<<< HEAD
    ret=open_files->next;


    /* check whether there were any to begin with */
    if (ret==open_files)
    {
        /* after this, the open_file pointer should never change */
        ret=NULL;
=======
    ret = open_files->next;


    /* check whether there were any to begin with */
    if (ret == open_files)
    {
        /* after this, the open_file pointer should never change */
        ret = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    else
    {
        gmx_fio_lock(open_files->next);
    }
    gmx_fio_unlock(open_files);


    return ret;
}

<<<<<<< HEAD
/* get the next open file, or NULL if there is none. 
=======
/* get the next open file, or NULL if there is none.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
   Unlocks the previous fio and locks the next one. */
static t_fileio *gmx_fio_get_next(t_fileio *fio)
{
    t_fileio *ret;

<<<<<<< HEAD
    ret=fio->next;
    /* check if that was the last one */
    if (fio->next==open_files)
    {
        ret=NULL;
#ifdef GMX_THREADS
=======
    ret = fio->next;
    /* check if that was the last one */
    if (fio->next == open_files)
    {
        ret = NULL;
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        tMPI_Thread_mutex_unlock(&open_file_mutex);
#endif
    }
    else
    {
        gmx_fio_lock(ret);
    }
    gmx_fio_unlock(fio);

    return ret;
}

/* Stop looping through the open_files.  Unlocks the global lock. */
static void gmx_fio_stop_getting_next(t_fileio *fio)
{
    gmx_fio_unlock(fio);
<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&open_file_mutex);
#endif
}




/*****************************************************************
 *
 *                     EXPORTED SECTION
 *
 *****************************************************************/
t_fileio *gmx_fio_open(const char *fn, const char *mode)
{
    t_fileio *fio = NULL;
<<<<<<< HEAD
    int i;
    char newmode[5];
    gmx_bool bRead, bReadWrite;
    int xdrid;
=======
    int       i;
    char      newmode[5];
    gmx_bool  bRead, bReadWrite;
    int       xdrid;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (fn2ftp(fn) == efTPA)
    {
        strcpy(newmode, mode);
    }
    else
    {
        /* sanitize the mode string */
        if (strncmp(mode, "r+", 2) == 0)
        {
            strcpy(newmode, "r+");
        }
        else if (mode[0] == 'r')
        {
            strcpy(newmode, "r");
        }
        else if (strncmp(mode, "w+", 2) == 0)
        {
            strcpy(newmode, "w+");
        }
        else if (mode[0] == 'w')
        {
            strcpy(newmode, "w");
        }
        else if (strncmp(mode, "a+", 2) == 0)
        {
            strcpy(newmode, "a+");
        }
        else if (mode[0] == 'a')
        {
            strcpy(newmode, "a");
        }
        else
        {
<<<<<<< HEAD
            gmx_fatal(FARGS, "DEATH HORROR in gmx_fio_open, mode is '%s'",mode);
=======
            gmx_fatal(FARGS, "DEATH HORROR in gmx_fio_open, mode is '%s'", mode);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }

    /* Check if it should be opened as a binary file */
<<<<<<< HEAD
    if (strncmp(ftp2ftype(fn2ftp(fn)),"ASCII",5))
    {
        /* Not ascii, add b to file mode */
        if ((strchr(newmode,'b')==NULL) && (strchr(newmode,'B')==NULL))
        {
            strcat(newmode,"b");
=======
    if (strncmp(ftp2ftype(fn2ftp(fn)), "ASCII", 5))
    {
        /* Not ascii, add b to file mode */
        if ((strchr(newmode, 'b') == NULL) && (strchr(newmode, 'B') == NULL))
        {
            strcat(newmode, "b");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }

    snew(fio, 1);
<<<<<<< HEAD
#ifdef GMX_THREADS
    tMPI_Lock_init(&(fio->mtx));
#endif
    bRead = (newmode[0]=='r' && newmode[1]!='+');
    bReadWrite = (newmode[1]=='+');
    fio->fp = NULL;
    fio->xdr = NULL;
    if (fn)
    {
        fio->iFTP = fn2ftp(fn);
        fio->fn = strdup(fn);
        fio->bStdio = FALSE;

        /* If this file type is in the list of XDR files, open it like that */
        if (in_ftpset(fio->iFTP,asize(ftpXDR),ftpXDR))
=======
#ifdef GMX_THREAD_MPI
    tMPI_Lock_init(&(fio->mtx));
#endif
    bRead      = (newmode[0] == 'r' && newmode[1] != '+');
    bReadWrite = (newmode[1] == '+');
    fio->fp    = NULL;
    fio->xdr   = NULL;
    if (fn)
    {
        fio->iFTP   = fn2ftp(fn);
        fio->fn     = strdup(fn);
        fio->bStdio = FALSE;

        /* If this file type is in the list of XDR files, open it like that */
        if (in_ftpset(fio->iFTP, asize(ftpXDR), ftpXDR))
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            /* First check whether we have to make a backup,
             * only for writing, not for read or append.
             */
<<<<<<< HEAD
            if (newmode[0]=='w')
=======
            if (newmode[0] == 'w')
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
#ifndef GMX_FAHCORE
                /* only make backups for normal gromacs */
                make_backup(fn);
#endif
            }
<<<<<<< HEAD
            else 
=======
            else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                /* Check whether file exists */
                if (!gmx_fexist(fn))
                {
                    gmx_open(fn);
                }
            }
            /* Open the file */
<<<<<<< HEAD
            fio->fp = ffopen(fn,newmode);

            /* determine the XDR direction */
            if (newmode[0] == 'w' || newmode[0]=='a')
            {
                fio->xdrmode=XDR_ENCODE;
            }
            else
            {
                fio->xdrmode=XDR_DECODE;
            }

            snew(fio->xdr,1);
=======
            fio->fp = ffopen(fn, newmode);

            /* determine the XDR direction */
            if (newmode[0] == 'w' || newmode[0] == 'a')
            {
                fio->xdrmode = XDR_ENCODE;
            }
            else
            {
                fio->xdrmode = XDR_DECODE;
            }

            snew(fio->xdr, 1);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            xdrstdio_create(fio->xdr, fio->fp, fio->xdrmode);
        }
        else
        {
            /* If it is not, open it as a regular file */
<<<<<<< HEAD
            fio->fp = ffopen(fn,newmode);
=======
            fio->fp = ffopen(fn, newmode);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }

        /* for appending seek to end of file to make sure ftell gives correct position
         * important for checkpointing */
<<<<<<< HEAD
        if (newmode[0]=='a')
=======
        if (newmode[0] == 'a')
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            gmx_fseek(fio->fp, 0, SEEK_END);
        }
    }
    else
    {
        /* Use stdin/stdout for I/O */
        fio->iFTP   = efTPA;
        fio->fp     = bRead ? stdin : stdout;
        fio->fn     = strdup("STDIO");
        fio->bStdio = TRUE;
    }
<<<<<<< HEAD
    fio->bRead  = bRead;
    fio->bReadWrite = bReadWrite;
    fio->bDouble= (sizeof(real) == sizeof(double));
    fio->bDebug = FALSE;
    fio->bOpen  = TRUE;
=======
    fio->bRead             = bRead;
    fio->bReadWrite        = bReadWrite;
    fio->bDouble           = (sizeof(real) == sizeof(double));
    fio->bDebug            = FALSE;
    fio->bOpen             = TRUE;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    fio->bLargerThan_off_t = FALSE;

    /* set the reader/writer functions */
    gmx_fio_set_iotype(fio);

    /* and now insert this file into the list of open files. */
    gmx_fio_insert(fio);
    return fio;
}

static int gmx_fio_close_locked(t_fileio *fio)
{
    int rc = 0;

    if (!fio->bOpen)
    {
<<<<<<< HEAD
        gmx_fatal(FARGS,"File %s closed twice!\n", fio->fn);
=======
        gmx_fatal(FARGS, "File %s closed twice!\n", fio->fn);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    if (in_ftpset(fio->iFTP, asize(ftpXDR), ftpXDR))
    {
        xdr_destroy(fio->xdr);
        sfree(fio->xdr);
    }

    /* Don't close stdin and stdout! */
<<<<<<< HEAD
    if (!fio->bStdio && fio->fp!=NULL)
        rc = ffclose(fio->fp); /* fclose returns 0 if happy */

=======
    if (!fio->bStdio && fio->fp != NULL)
    {
        rc = ffclose(fio->fp); /* fclose returns 0 if happy */

    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    fio->bOpen = FALSE;

    return rc;
}

int gmx_fio_close(t_fileio *fio)
{
    int rc = 0;

<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* first lock the big open_files mutex. */
    /* We don't want two processes operating on the list at the same time */
    tMPI_Thread_mutex_lock(&open_file_mutex);
#endif

    gmx_fio_lock(fio);
    /* first remove it from the list */
    gmx_fio_remove(fio);
<<<<<<< HEAD
    rc=gmx_fio_close_locked(fio);
=======
    rc = gmx_fio_close_locked(fio);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    sfree(fio);

<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&open_file_mutex);
#endif

    return rc;
}

/* close only fp but keep FIO entry. */
int gmx_fio_fp_close(t_fileio *fio)
{
<<<<<<< HEAD
    int rc=0;
    gmx_fio_lock(fio);
    if (!in_ftpset(fio->iFTP,asize(ftpXDR),ftpXDR) && !fio->bStdio)
    {
        rc = ffclose(fio->fp); /* fclose returns 0 if happy */
        fio->fp = NULL; 
=======
    int rc = 0;
    gmx_fio_lock(fio);
    if (!in_ftpset(fio->iFTP, asize(ftpXDR), ftpXDR) && !fio->bStdio)
    {
        rc      = ffclose(fio->fp); /* fclose returns 0 if happy */
        fio->fp = NULL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    gmx_fio_unlock(fio);

    return rc;
}

FILE * gmx_fio_fopen(const char *fn, const char *mode)
{
<<<<<<< HEAD
    FILE *fp, *ret;
=======
    FILE     *fp, *ret;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_fileio *fio;

    fio = gmx_fio_open(fn, mode);
    gmx_fio_lock(fio);
    ret = fio->fp;
    gmx_fio_unlock(fio);

    return ret;
}

int gmx_fio_fclose(FILE *fp)
{
    t_fileio *cur;
<<<<<<< HEAD
    t_fileio *found=NULL;
    int rc=-1;

    cur=gmx_fio_get_first();
    while(cur)
    {
        if (cur->fp == fp)
        {
            rc=gmx_fio_close_locked(cur);
=======
    t_fileio *found = NULL;
    int       rc    = -1;

    cur = gmx_fio_get_first();
    while (cur)
    {
        if (cur->fp == fp)
        {
            rc = gmx_fio_close_locked(cur);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            gmx_fio_remove(cur);
            gmx_fio_stop_getting_next(cur);
            break;
        }
<<<<<<< HEAD
        cur=gmx_fio_get_next(cur);
=======
        cur = gmx_fio_get_next(cur);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return rc;
}

/* internal variant of get_file_md5 that operates on a locked file */
<<<<<<< HEAD
static int gmx_fio_int_get_file_md5(t_fileio *fio, gmx_off_t offset, 
                                    unsigned char digest[])
{
    /*1MB: large size important to catch almost identical files */
#define CPT_CHK_LEN  1048576 
    md5_state_t state;
    unsigned char buf[CPT_CHK_LEN];
    gmx_off_t read_len;
    gmx_off_t seek_offset;
    int ret = -1;
=======
static int gmx_fio_int_get_file_md5(t_fileio *fio, gmx_off_t offset,
                                    unsigned char digest[])
{
    /*1MB: large size important to catch almost identical files */
#define CPT_CHK_LEN  1048576
    md5_state_t   state;
    unsigned char buf[CPT_CHK_LEN];
    gmx_off_t     read_len;
    gmx_off_t     seek_offset;
    int           ret = -1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    seek_offset = offset - CPT_CHK_LEN;
    if (seek_offset < 0)
    {
        seek_offset = 0;
    }
    read_len = offset - seek_offset;


    if (fio->fp && fio->bReadWrite)
    {
<<<<<<< HEAD
        ret=gmx_fseek(fio->fp, seek_offset, SEEK_SET);
=======
        ret = gmx_fseek(fio->fp, seek_offset, SEEK_SET);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (ret)
        {
            gmx_fseek(fio->fp, 0, SEEK_END);
        }
    }
    if (ret) /*either no fp, not readwrite, or fseek not successful */
    {
        return -1;
    }

    /* the read puts the file position back to offset */
    if ((gmx_off_t)fread(buf, 1, read_len, fio->fp) != read_len)
    {
        /* not fatal: md5sum check to prevent overwriting files
         * works (less safe) without
         * */
        if (ferror(fio->fp))
        {
            fprintf(stderr, "\nTrying to get md5sum: %s: %s\n", fio->fn,
                    strerror(errno));
        }
        else if (feof(fio->fp))
        {
            /*
             * For long runs that checkpoint frequently but write e.g. logs
             * infrequently we don't want to issue lots of warnings before we
             * have written anything to the log.
             */
<<<<<<< HEAD
            if(0)
=======
            if (0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                fprintf(stderr, "\nTrying to get md5sum: EOF: %s\n", fio->fn);
            }
        }
        else
        {
            fprintf(
<<<<<<< HEAD
                stderr,
                "\nTrying to get md5sum: Unknown reason for short read: %s\n",
                fio->fn);
=======
                    stderr,
                    "\nTrying to get md5sum: Unknown reason for short read: %s\n",
                    fio->fn);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }

        gmx_fseek(fio->fp, 0, SEEK_END);

        ret = -1;
    }
<<<<<<< HEAD
    gmx_fseek(fio->fp, 0, SEEK_END); /*is already at end, but under windows 
                                       it gives problems otherwise*/
=======
    gmx_fseek(fio->fp, 0, SEEK_END); /*is already at end, but under windows
                                        it gives problems otherwise*/
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (debug)
    {
        fprintf(debug, "chksum %s readlen %ld\n", fio->fn, (long int)read_len);
    }

    if (!ret)
    {
        md5_init(&state);
        md5_append(&state, buf, read_len);
        md5_finish(&state, digest);
        return read_len;
    }
    else
    {
        return ret;
    }
}


/*
 * fio: file to compute md5 for
 * offset: starting pointer of region to use for md5
<<<<<<< HEAD
 * digest: return array of md5 sum 
 */
int gmx_fio_get_file_md5(t_fileio *fio, gmx_off_t offset, 
=======
 * digest: return array of md5 sum
 */
int gmx_fio_get_file_md5(t_fileio *fio, gmx_off_t offset,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                         unsigned char digest[])
{
    int ret;

    gmx_fio_lock(fio);
<<<<<<< HEAD
    ret=gmx_fio_int_get_file_md5(fio, offset, digest);
=======
    ret = gmx_fio_int_get_file_md5(fio, offset, digest);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    return ret;
}

/* The fio_mutex should ALWAYS be locked when this function is called */
static int gmx_fio_int_get_file_position(t_fileio *fio, gmx_off_t *offset)
{
    char buf[STRLEN];

    /* Flush the file, so we are sure it is written */
    if (gmx_fio_int_flush(fio))
    {
        char buf[STRLEN];
        sprintf(
<<<<<<< HEAD
            buf,
            "Cannot write file '%s'; maybe you are out of disk space or quota?",
            fio->fn);
        gmx_file(buf);
    }

    /* We cannot count on XDR being able to write 64-bit integers, 
     so separate into high/low 32-bit values.
     In case the filesystem has 128-bit offsets we only care 
     about the first 64 bits - we'll have to fix
     this when exabyte-size output files are common...
     */
    *offset=gmx_ftell(fio->fp);
=======
                buf,
                "Cannot write file '%s'; maybe you are out of disk space?",
                fio->fn);
        gmx_file(buf);
    }

    /* We cannot count on XDR being able to write 64-bit integers,
       so separate into high/low 32-bit values.
       In case the filesystem has 128-bit offsets we only care
       about the first 64 bits - we'll have to fix
       this when exabyte-size output files are common...
     */
    *offset = gmx_ftell(fio->fp);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return 0;
}

int gmx_fio_check_file_position(t_fileio *fio)
{
    /* If gmx_off_t is 4 bytes we can not store file offset > 2 GB.
     * If we do not have ftello, we will play it safe.
     */
#if (SIZEOF_GMX_OFF_T == 4 || !defined HAVE_FSEEKO)
    gmx_off_t offset;

    gmx_fio_lock(fio);
<<<<<<< HEAD
    gmx_fio_int_get_file_position(fio,&offset);
=======
    gmx_fio_int_get_file_position(fio, &offset);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* We have a 4 byte offset,
     * make sure that we will detect out of range for all possible cases.
     */
    if (offset < 0 || offset > 2147483647)
    {
        fio->bLargerThan_off_t = TRUE;
    }
    gmx_fio_unlock(fio);
#endif

    return 0;
}

int gmx_fio_get_output_file_positions(gmx_file_position_t **p_outputfiles,
<<<<<<< HEAD
                                      int *p_nfiles)
{
    int i, nfiles, rc, nalloc;
    int pos_hi, pos_lo;
    long pos;
    gmx_file_position_t * outputfiles;
    char buf[STRLEN];
    t_fileio *cur;
=======
                                      int                  *p_nfiles)
{
    int                   i, nfiles, rc, nalloc;
    int                   pos_hi, pos_lo;
    long                  pos;
    gmx_file_position_t * outputfiles;
    char                  buf[STRLEN];
    t_fileio             *cur;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    nfiles = 0;

    /* pre-allocate 100 files */
    nalloc = 100;
<<<<<<< HEAD
    snew(outputfiles,nalloc);

    cur=gmx_fio_get_first();
    while(cur)
    {
        /* Skip the checkpoint files themselves, since they could be open when 
           we call this routine... */
        /* also skip debug files (shoud be the only iFTP==efNR) */
        if (cur->bOpen && 
            !cur->bRead && 
            !cur->bStdio && 
            cur->iFTP != efCPT && 
=======
    snew(outputfiles, nalloc);

    cur = gmx_fio_get_first();
    while (cur)
    {
        /* Skip the checkpoint files themselves, since they could be open when
           we call this routine... */
        /* also skip debug files (shoud be the only iFTP==efNR) */
        if (cur->bOpen &&
            !cur->bRead &&
            !cur->bStdio &&
            cur->iFTP != efCPT &&
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            cur->iFTP != efNR)
        {
            int ret;
            /* This is an output file currently open for writing, add it */
            if (nfiles == nalloc)
            {
                nalloc += 100;
<<<<<<< HEAD
                srenew(outputfiles,nalloc);
=======
                srenew(outputfiles, nalloc);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }

            strncpy(outputfiles[nfiles].filename, cur->fn, STRLEN - 1);

            /* Get the file position */
            if (cur->bLargerThan_off_t)
            {
                /* -1 signals out of range */
<<<<<<< HEAD
                outputfiles[nfiles].offset = -1;
=======
                outputfiles[nfiles].offset      = -1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                outputfiles[nfiles].chksum_size = -1;
            }
            else
            {
                gmx_fio_int_get_file_position(cur, &outputfiles[nfiles].offset);
#ifndef GMX_FAHCORE
                outputfiles[nfiles].chksum_size
<<<<<<< HEAD
                    = gmx_fio_int_get_file_md5(cur, 
=======
                    = gmx_fio_int_get_file_md5(cur,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                                               outputfiles[nfiles].offset,
                                               outputfiles[nfiles].chksum);
#endif
            }

            nfiles++;
        }

<<<<<<< HEAD
        cur=gmx_fio_get_next(cur);
    }
    *p_nfiles = nfiles;
=======
        cur = gmx_fio_get_next(cur);
    }
    *p_nfiles      = nfiles;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    *p_outputfiles = outputfiles;

    return 0;
}


void gmx_fio_checktype(t_fileio *fio)
{
    if (in_ftpset(fio->iFTP, asize(ftpXDR), ftpXDR))
    {
        return;
    }
    else if (in_ftpset(fio->iFTP, asize(ftpASC), ftpASC))
    {
        return;
    }
    else if (in_ftpset(fio->iFTP, asize(ftpBIN), ftpBIN))
    {
        return;
    }
<<<<<<< HEAD
#ifdef HAVE_XMl
    else if (in_ftpset(fio->iFTP,asize(ftpXML),ftpXML))
    {
        return;
    }
#endif
    else
        gmx_fatal(FARGS, "Can not read/write topologies to file type %s",
                  ftp2ext(fio->iFTP));
=======
    else
    {
        gmx_fatal(FARGS, "Can not read/write topologies to file type %s",
                  ftp2ext(fio->iFTP));
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

}


void gmx_fio_setprecision(t_fileio *fio, gmx_bool bDouble)
{
    gmx_fio_lock(fio);
    fio->bDouble = bDouble;
    gmx_fio_unlock(fio);
}

gmx_bool gmx_fio_getdebug(t_fileio *fio)
{
    gmx_bool ret;

    gmx_fio_lock(fio);
    ret = fio->bDebug;
    gmx_fio_unlock(fio);

    return ret;
}

void gmx_fio_setdebug(t_fileio *fio, gmx_bool bDebug)
{
    gmx_fio_lock(fio);
    fio->bDebug = bDebug;
    gmx_fio_unlock(fio);
}

char *gmx_fio_getname(t_fileio *fio)
{
    char *ret;
    gmx_fio_lock(fio);
    ret = fio->fn;
    gmx_fio_unlock(fio);

    return ret;
}

int gmx_fio_getftp(t_fileio* fio)
{
    int ret;

    gmx_fio_lock(fio);
    ret = fio->iFTP;
    gmx_fio_unlock(fio);

    return ret;
}

void gmx_fio_rewind(t_fileio* fio)
{
    gmx_fio_lock(fio);

    if (fio->xdr)
    {
        xdr_destroy(fio->xdr);
        frewind(fio->fp);
        xdrstdio_create(fio->xdr, fio->fp, fio->xdrmode);
    }
    else
    {
        frewind(fio->fp);
    }
    gmx_fio_unlock(fio);
}


int gmx_fio_flush(t_fileio* fio)
{
    int ret;

    gmx_fio_lock(fio);
<<<<<<< HEAD
    ret=gmx_fio_int_flush(fio);
=======
    ret = gmx_fio_int_flush(fio);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    return ret;
}



static int gmx_fio_int_fsync(t_fileio *fio)
{
<<<<<<< HEAD
    int rc = 0;
    int filen=-1;
=======
    int rc    = 0;
    int filen = -1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    if (fio->fp)
    {
<<<<<<< HEAD
        rc=gmx_fsync(fio->fp);
    }
    else if (fio->xdr) /* this should normally not happen */
    {
        rc=gmx_fsync((FILE*) fio->xdr->x_private);
                                    /* ^ is this actually OK? */
=======
        rc = gmx_fsync(fio->fp);
    }
    else if (fio->xdr) /* this should normally not happen */
    {
        rc = gmx_fsync((FILE*) fio->xdr->x_private);
        /* ^ is this actually OK? */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    return rc;
}


int gmx_fio_fsync(t_fileio *fio)
{
    int rc;

    gmx_fio_lock(fio);
<<<<<<< HEAD
    rc=gmx_fio_int_fsync(fio);
=======
    rc = gmx_fio_int_fsync(fio);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    return rc;
}



t_fileio *gmx_fio_all_output_fsync(void)
{
<<<<<<< HEAD
    t_fileio *ret=NULL;
    t_fileio *cur;

    cur=gmx_fio_get_first();
    while(cur)
    {
        /* skip debug files (shoud be the only iFTP==efNR) */
        if (cur->bOpen && 
            !cur->bRead && 
            !cur->bStdio && 
            cur->iFTP != efNR)
        {
            /* if any of them fails, return failure code */
            int rc=gmx_fio_int_fsync(cur);
            if (rc != 0 && !ret) 
            {
                ret=cur;
            }
        }
        cur=gmx_fio_get_next(cur);
=======
    t_fileio *ret = NULL;
    t_fileio *cur;

    cur = gmx_fio_get_first();
    while (cur)
    {
        /* skip debug files (shoud be the only iFTP==efNR) */
        if (cur->bOpen &&
            !cur->bRead &&
            !cur->bStdio &&
            cur->iFTP != efNR)
        {
            /* if any of them fails, return failure code */
            int rc = gmx_fio_int_fsync(cur);
            if (rc != 0 && !ret)
            {
                ret = cur;
            }
        }
        cur = gmx_fio_get_next(cur);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    /* in addition, we force these to be written out too, if they're being
       redirected. We don't check for errors because errors most likely mean
       that they're not redirected. */
    fflush(stdout);
    fflush(stderr);
#if (defined(HAVE_FSYNC))
    /* again, fahcore defines HAVE_FSYNC and fsync() */
    fsync(STDOUT_FILENO);
    fsync(STDERR_FILENO);
#endif

    return ret;
}


gmx_off_t gmx_fio_ftell(t_fileio* fio)
{
    gmx_off_t ret = 0;

    gmx_fio_lock(fio);
    if (fio->fp)
<<<<<<< HEAD
        ret = gmx_ftell(fio->fp);
=======
    {
        ret = gmx_ftell(fio->fp);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);
    return ret;
}

int gmx_fio_seek(t_fileio* fio, gmx_off_t fpos)
{
    int rc;

    gmx_fio_lock(fio);
    if (fio->fp)
    {
        rc = gmx_fseek(fio->fp, fpos, SEEK_SET);
    }
    else
    {
        gmx_file(fio->fn);
        rc = -1;
    }
    gmx_fio_unlock(fio);
    return rc;
}

FILE *gmx_fio_getfp(t_fileio *fio)
{
    FILE *ret = NULL;

    gmx_fio_lock(fio);
    if (fio->fp)
<<<<<<< HEAD
        ret = fio->fp;
=======
    {
        ret = fio->fp;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);
    return ret;
}

XDR *gmx_fio_getxdr(t_fileio* fio)
{
    XDR *ret = NULL;

    gmx_fio_lock(fio);
    if (fio->xdr)
<<<<<<< HEAD
        ret = fio->xdr;
=======
    {
        ret = fio->xdr;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    return ret;
}

gmx_bool gmx_fio_getread(t_fileio* fio)
{
    gmx_bool ret;

    gmx_fio_lock(fio);
    ret = fio->bRead;
    gmx_fio_unlock(fio);

    return ret;
}

int xtc_seek_frame(t_fileio *fio, int frame, int natoms)
{
    int ret;

    gmx_fio_lock(fio);
<<<<<<< HEAD
    ret=xdr_xtc_seek_frame(frame, fio->fp, fio->xdr, natoms);
=======
    ret = xdr_xtc_seek_frame(frame, fio->fp, fio->xdr, natoms);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    return ret;
}

<<<<<<< HEAD
int xtc_seek_time(t_fileio *fio, real time, int natoms)
=======
int xtc_seek_time(t_fileio *fio, real time, int natoms, gmx_bool bSeekForwardOnly)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    int ret;

    gmx_fio_lock(fio);
<<<<<<< HEAD
    ret=xdr_xtc_seek_time(time, fio->fp, fio->xdr, natoms);
=======
    ret = xdr_xtc_seek_time(time, fio->fp, fio->xdr, natoms, bSeekForwardOnly);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_fio_unlock(fio);

    return ret;
}
<<<<<<< HEAD


=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
