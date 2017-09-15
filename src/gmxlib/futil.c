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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef HAVE_DIRENT_H
/* POSIX */
#include <dirent.h>
#endif


<<<<<<< HEAD
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
=======
#ifdef GMX_NATIVE_WINDOWS
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include <direct.h>
#include <io.h>
#endif

#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "network.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "statutil.h"


<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "thread_mpi.h"
#endif

/* Windows file stuff, only necessary for visual studio */
#ifdef _MSC_VER
#include "windows.h"
#endif

<<<<<<< HEAD
/* we keep a linked list of all files opened through pipes (i.e. 
   compressed or .gzipped files. This way we can distinguish between them
   without having to change the semantics of reading from/writing to files) 
   */
typedef struct t_pstack {
    FILE   *fp;
    struct t_pstack *prev;
} t_pstack;

static t_pstack *pstack=NULL;
static gmx_bool     bUnbuffered=FALSE;

#ifdef GMX_THREADS
/* this linked list is an intrinsically globally shared object, so we have
   to protect it with mutexes */
static tMPI_Thread_mutex_t pstack_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
=======
/* we keep a linked list of all files opened through pipes (i.e.
   compressed or .gzipped files. This way we can distinguish between them
   without having to change the semantics of reading from/writing to files)
 */
typedef struct t_pstack {
    FILE            *fp;
    struct t_pstack *prev;
} t_pstack;

static t_pstack    *pstack      = NULL;
static gmx_bool     bUnbuffered = FALSE;

#ifdef GMX_THREAD_MPI
/* this linked list is an intrinsically globally shared object, so we have
   to protect it with mutexes */
static tMPI_Thread_mutex_t pstack_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

void no_buffers(void)
{
<<<<<<< HEAD
    bUnbuffered=TRUE;
=======
    bUnbuffered = TRUE;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

void push_ps(FILE *fp)
{
    t_pstack *ps;

<<<<<<< HEAD
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    snew(ps,1);
    ps->fp   = fp;
    ps->prev = pstack;
    pstack   = ps;
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    snew(ps, 1);
    ps->fp   = fp;
    ps->prev = pstack;
    pstack   = ps;
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
}

#ifdef GMX_FAHCORE
/* don't use pipes!*/
#define popen fah_fopen
#define pclose fah_fclose
#define SKIP_FFOPS 1
#else
#ifdef ffclose
#undef ffclose
#endif
#endif

#ifndef GMX_FAHCORE
#ifndef HAVE_PIPES
<<<<<<< HEAD
static FILE *popen(const char *nm,const char *mode)
=======
static FILE *popen(const char *nm, const char *mode)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    gmx_impl("Sorry no pipes...");

    return NULL;
}

static int pclose(FILE *fp)
{
    gmx_impl("Sorry no pipes...");

    return 0;
}
#endif
#endif

int ffclose(FILE *fp)
{
#ifdef SKIP_FFOPS
    return fclose(fp);
#else
<<<<<<< HEAD
    t_pstack *ps,*tmp;
    int ret=0;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps=pstack;
    if (ps == NULL) {
        if (fp != NULL) 
            ret = fclose(fp);
    }
    else if (ps->fp == fp) {
        if (fp != NULL)
            ret = pclose(fp);
        pstack=pstack->prev;
        sfree(ps);
    }
    else {
        while ((ps->prev != NULL) && (ps->prev->fp != fp))
            ps=ps->prev;
        if (ps->prev->fp == fp) {
            if (ps->prev->fp != NULL)
                ret = pclose(ps->prev->fp);
            tmp=ps->prev;
            ps->prev=ps->prev->prev;
            sfree(tmp);
        }
        else {
            if (fp != NULL)
                ret = fclose(fp);
        }
    }
#ifdef GMX_THREADS
=======
    t_pstack *ps, *tmp;
    int       ret = 0;
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps = pstack;
    if (ps == NULL)
    {
        if (fp != NULL)
        {
            ret = fclose(fp);
        }
    }
    else if (ps->fp == fp)
    {
        if (fp != NULL)
        {
            ret = pclose(fp);
        }
        pstack = pstack->prev;
        sfree(ps);
    }
    else
    {
        while ((ps->prev != NULL) && (ps->prev->fp != fp))
        {
            ps = ps->prev;
        }
        if ((ps->prev != NULL) && ps->prev->fp == fp)
        {
            if (ps->prev->fp != NULL)
            {
                ret = pclose(ps->prev->fp);
            }
            tmp      = ps->prev;
            ps->prev = ps->prev->prev;
            sfree(tmp);
        }
        else
        {
            if (fp != NULL)
            {
                ret = fclose(fp);
            }
        }
    }
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
    return ret;
#endif
}


#ifdef rewind
#undef rewind
#endif

void frewind(FILE *fp)
{
    t_pstack *ps;
<<<<<<< HEAD
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps=pstack;
    while (ps != NULL) {
        if (ps->fp == fp) {
            fprintf(stderr,"Cannot rewind compressed file!\n");
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps = pstack;
    while (ps != NULL)
    {
        if (ps->fp == fp)
        {
            fprintf(stderr, "Cannot rewind compressed file!\n");
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
            return;
        }
<<<<<<< HEAD
        ps=ps->prev;
    }
    rewind(fp);
#ifdef GMX_THREADS
=======
        ps = ps->prev;
    }
    rewind(fp);
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
}

int gmx_fseek(FILE *stream, gmx_off_t offset, int whence)
{
#ifdef HAVE_FSEEKO
    return fseeko(stream, offset, whence);
#else
#ifdef HAVE__FSEEKI64
    return _fseeki64(stream, offset, whence);
#else
    return fseek(stream, offset, whence);
#endif
#endif
}

gmx_off_t gmx_ftell(FILE *stream)
{
#ifdef HAVE_FSEEKO
    return ftello(stream);
#else
<<<<<<< HEAD
#ifdef HAVE__FSEEKI64 
=======
#ifdef HAVE__FSEEKI64
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return _ftelli64(stream);
#else
    return ftell(stream);
#endif
#endif
}


gmx_bool is_pipe(FILE *fp)
{
    t_pstack *ps;
<<<<<<< HEAD
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps=pstack;
    while (ps != NULL) {
        if (ps->fp == fp) {
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps = pstack;
    while (ps != NULL)
    {
        if (ps->fp == fp)
        {
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
            return TRUE;
        }
<<<<<<< HEAD
        ps=ps->prev;
    }
#ifdef GMX_THREADS
=======
        ps = ps->prev;
    }
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
    return FALSE;
}


<<<<<<< HEAD
static FILE *uncompress(const char *fn,const char *mode)
{
    FILE *fp;
    char buf[256];

    sprintf(buf,"uncompress -c < %s",fn);
    fprintf(stderr,"Going to execute '%s'\n",buf);
    if ((fp=popen(buf,mode)) == NULL)
        gmx_open(fn);
=======
static FILE *uncompress(const char *fn, const char *mode)
{
    FILE *fp;
    char  buf[256];

    sprintf(buf, "uncompress -c < %s", fn);
    fprintf(stderr, "Going to execute '%s'\n", buf);
    if ((fp = popen(buf, mode)) == NULL)
    {
        gmx_open(fn);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    push_ps(fp);

    return fp;
}

<<<<<<< HEAD
static FILE *gunzip(const char *fn,const char *mode)
{
    FILE *fp;
    char buf[256];

    sprintf(buf,"gunzip -c < %s",fn);
    fprintf(stderr,"Going to execute '%s'\n",buf);
    if ((fp=popen(buf,mode)) == NULL)
        gmx_open(fn);
=======
static FILE *gunzip(const char *fn, const char *mode)
{
    FILE *fp;
    char  buf[256];

    sprintf(buf, "gunzip -c < %s", fn);
    fprintf(stderr, "Going to execute '%s'\n", buf);
    if ((fp = popen(buf, mode)) == NULL)
    {
        gmx_open(fn);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    push_ps(fp);

    return fp;
}

gmx_bool gmx_fexist(const char *fname)
{
    FILE *test;

    if (fname == NULL)
<<<<<<< HEAD
        return FALSE;
    test=fopen(fname,"r");
    if (test == NULL) {
        /*Windows doesn't allow fopen of directory - so we need to check this seperately */
        #if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__) 
            DWORD attr = GetFileAttributes(fname);
            return (attr != INVALID_FILE_ATTRIBUTES) && (attr & FILE_ATTRIBUTE_DIRECTORY);
        #else 
            return FALSE;
        #endif
    } else {
=======
    {
        return FALSE;
    }
    test = fopen(fname, "r");
    if (test == NULL)
    {
        /*Windows doesn't allow fopen of directory - so we need to check this seperately */
        #ifdef GMX_NATIVE_WINDOWS
        DWORD attr = GetFileAttributes(fname);
        return (attr != INVALID_FILE_ATTRIBUTES) && (attr & FILE_ATTRIBUTE_DIRECTORY);
        #else
        return FALSE;
        #endif
    }
    else
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fclose(test);
        return TRUE;
    }
}

<<<<<<< HEAD

gmx_bool gmx_fexist_master(const char *fname, t_commrec *cr)
{
  gmx_bool bExist;
  
  if (SIMMASTER(cr)) 
  {
      bExist = gmx_fexist(fname);
  }
  if (PAR(cr)) 
  {
      gmx_bcast(sizeof(bExist),&bExist,cr);
  }
  return bExist;
=======
static gmx_bool gmx_is_file(const char *fname)
{
    FILE *test;

    if (fname == NULL)
    {
        return FALSE;
    }
    test = fopen(fname, "r");
    if (test == NULL)
    {
        return FALSE;
    }
    else
    {
        fclose(test);
        /*Windows doesn't allow fopen of directory - so we don't need to check this seperately */
        #if (!((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__))
        {
            int         status;
            struct stat st_buf;
            status = stat (fname, &st_buf);
            if (status != 0 || !S_ISREG(st_buf.st_mode))
            {
                return FALSE;
            }
        }
        #endif
        return TRUE;
    }
}


gmx_bool gmx_fexist_master(const char *fname, t_commrec *cr)
{
    gmx_bool bExist;

    if (SIMMASTER(cr))
    {
        bExist = gmx_fexist(fname);
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(bExist), &bExist, cr);
    }
    return bExist;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

gmx_bool gmx_eof(FILE *fp)
{
<<<<<<< HEAD
    char data[4];
    gmx_bool beof;

    if (is_pipe(fp))
        return feof(fp);
    else {
        if ((beof=fread(data,1,1,fp))==1)
            gmx_fseek(fp,-1,SEEK_CUR);
=======
    char     data[4];
    gmx_bool beof;

    if (is_pipe(fp))
    {
        return feof(fp);
    }
    else
    {
        if ((beof = fread(data, 1, 1, fp)) == 1)
        {
            gmx_fseek(fp, -1, SEEK_CUR);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        return !beof;
    }
}

<<<<<<< HEAD
static char *backup_fn(const char *file,int count_max)
=======
static char *backup_fn(const char *file, int count_max)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    /* Use a reasonably low value for countmax; we might
     * generate 4-5 files in each round, and we dont
     * want to hit directory limits of 1024 or 2048 files.
     */
#define COUNTMAX 99
<<<<<<< HEAD
    int         i,count=1;
    char        *directory,*fn;
=======
    int          i, count = 1;
    char        *directory, *fn;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    char        *buf;

    if (count_max == -1)
    {
        count_max = COUNTMAX;
    }

    smalloc(buf, GMX_PATH_MAX);

<<<<<<< HEAD
    for(i=strlen(file)-1; ((i > 0) && (file[i] != DIR_SEPARATOR)); i--)
        ;
=======
    for (i = strlen(file)-1; ((i > 0) && (file[i] != DIR_SEPARATOR)); i--)
    {
        ;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    /* Must check whether i > 0, i.e. whether there is a directory
     * in the file name. In that case we overwrite the / sign with
     * a '\0' to end the directory string .
     */
<<<<<<< HEAD
    if (i > 0) {
=======
    if (i > 0)
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        directory    = gmx_strdup(file);
        directory[i] = '\0';
        fn           = gmx_strdup(file+i+1);
    }
<<<<<<< HEAD
    else {
        directory    = gmx_strdup(".");
        fn           = gmx_strdup(file);
    }
    do {
        sprintf(buf,"%s/#%s.%d#",directory,fn,count);
        count++;
    } while ((count <= count_max) && gmx_fexist(buf));

    /* Arbitrarily bail out */
    if (count > count_max) 
        gmx_fatal(FARGS,"Won't make more than %d backups of %s for you.\n"
                  "The env.var. GMX_MAXBACKUP controls this maximum, -1 disables backups.",
                  count_max,fn);
=======
    else
    {
        directory    = gmx_strdup(".");
        fn           = gmx_strdup(file);
    }
    do
    {
        sprintf(buf, "%s/#%s.%d#", directory, fn, count);
        count++;
    }
    while ((count <= count_max) && gmx_fexist(buf));

    /* Arbitrarily bail out */
    if (count > count_max)
    {
        gmx_fatal(FARGS, "Won't make more than %d backups of %s for you.\n"
                  "The env.var. GMX_MAXBACKUP controls this maximum, -1 disables backups.",
                  count_max, fn);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    sfree(directory);
    sfree(fn);

    return buf;
}

gmx_bool make_backup(const char * name)
{
    char * env;
<<<<<<< HEAD
    int  count_max;
=======
    int    count_max;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    char * backup;

#ifdef GMX_FAHCORE
    return FALSE; /* skip making backups */
#else

    if (gmx_fexist(name))
    {
        env = getenv("GMX_MAXBACKUP");
        if (env != NULL)
        {
            count_max = 0;
<<<<<<< HEAD
            sscanf(env,"%d",&count_max);
=======
            sscanf(env, "%d", &count_max);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (count_max == -1)
            {
                /* Do not make backups and possibly overwrite old files */
                return TRUE;
            }
        }
        else
        {
            /* Use the default maximum */
            count_max = -1;
        }
<<<<<<< HEAD
        backup = backup_fn(name,count_max);
        if(rename(name, backup) == 0) {
            fprintf(stderr, "\nBack Off! I just backed up %s to %s\n",
                    name, backup);
        } else {
=======
        backup = backup_fn(name, count_max);
        if (rename(name, backup) == 0)
        {
            fprintf(stderr, "\nBack Off! I just backed up %s to %s\n",
                    name, backup);
        }
        else
        {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            fprintf(stderr, "Sorry couldn't backup %s to %s\n", name, backup);
            return FALSE;
        }
        sfree(backup);
    }
    return TRUE;
#endif
}

<<<<<<< HEAD
FILE *ffopen(const char *file,const char *mode)
{
#ifdef SKIP_FFOPS
    return fopen(file,mode);
#else
    FILE *ff=NULL;
    char buf[256],*bf,*bufsize=0,*ptr;
    gmx_bool bRead;
    int  bs;

    if (mode[0]=='w') {
=======
FILE *ffopen(const char *file, const char *mode)
{
#ifdef SKIP_FFOPS
    return fopen(file, mode);
#else
    FILE    *ff = NULL;
    char     buf[256], *bf, *bufsize = 0, *ptr;
    gmx_bool bRead;
    int      bs;

    if (file == NULL)
    {
        return NULL;
    }

    if (mode[0] == 'w')
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        make_backup(file);
    }
    where();

<<<<<<< HEAD
    bRead= (mode[0]=='r'&&mode[1]!='+');
    strcpy(buf,file);
    if (!bRead || gmx_fexist(buf)) {
        if ((ff=fopen(buf,mode))==NULL)
            gmx_file(buf);
=======
    bRead = (mode[0] == 'r' && mode[1] != '+');
    strcpy(buf, file);
    if (!bRead || gmx_fexist(buf))
    {
        if ((ff = fopen(buf, mode)) == NULL)
        {
            gmx_file(buf);
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        where();
        /* Check whether we should be using buffering (default) or not
         * (for debugging)
         */
<<<<<<< HEAD
        if (bUnbuffered || ((bufsize=getenv("LOG_BUFS")) != NULL)) {
            /* Check whether to use completely unbuffered */
            if (bUnbuffered)
                bs = 0;
            else
                bs=strtol(bufsize, NULL, 10); 
            if (bs <= 0)
                setbuf(ff,NULL); 
            else {
                snew(ptr,bs+8);
                if (setvbuf(ff,ptr,_IOFBF,bs) != 0)
                    gmx_file("Buffering File");
=======
        if (bUnbuffered || ((bufsize = getenv("LOG_BUFS")) != NULL))
        {
            /* Check whether to use completely unbuffered */
            if (bUnbuffered)
            {
                bs = 0;
            }
            else
            {
                bs = strtol(bufsize, NULL, 10);
            }
            if (bs <= 0)
            {
                setbuf(ff, NULL);
            }
            else
            {
                snew(ptr, bs+8);
                if (setvbuf(ff, ptr, _IOFBF, bs) != 0)
                {
                    gmx_file("Buffering File");
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
        where();
    }
<<<<<<< HEAD
    else {
        sprintf(buf,"%s.Z",file);
        if (gmx_fexist(buf)) {
            ff=uncompress(buf,mode);
        }
        else {
            sprintf(buf,"%s.gz",file);
            if (gmx_fexist(buf)) {
                ff=gunzip(buf,mode);
            }
            else 
                gmx_file(file);
=======
    else
    {
        sprintf(buf, "%s.Z", file);
        if (gmx_fexist(buf))
        {
            ff = uncompress(buf, mode);
        }
        else
        {
            sprintf(buf, "%s.gz", file);
            if (gmx_fexist(buf))
            {
                ff = gunzip(buf, mode);
            }
            else
            {
                gmx_file(file);
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
    return ff;
#endif
}

/* Our own implementation of dirent-like functionality to scan directories. */
struct gmx_directory
{
#ifdef HAVE_DIRENT_H
    DIR  *               dirent_handle;
<<<<<<< HEAD
#elif (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64)
=======
#elif (defined GMX_NATIVE_WINDOWS)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    intptr_t             windows_handle;
    struct _finddata_t   finddata;
    int                  first;
#else
<<<<<<< HEAD
    int      dummy;
=======
    int                  dummy;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
};


int
<<<<<<< HEAD
gmx_directory_open(gmx_directory_t *p_gmxdir,const char *dirname)
{
    struct gmx_directory *  gmxdir;
    int                     rc;
    
    snew(gmxdir,1);
    
    *p_gmxdir = gmxdir;
    
#ifdef HAVE_DIRENT_H
    if( (gmxdir->dirent_handle = opendir(dirname)) != NULL)
    {
        rc = 0;
    }
    else 
=======
gmx_directory_open(gmx_directory_t *p_gmxdir, const char *dirname)
{
    struct gmx_directory *  gmxdir;
    int                     rc;

    snew(gmxdir, 1);

    *p_gmxdir = gmxdir;

#ifdef HAVE_DIRENT_H
    if ( (gmxdir->dirent_handle = opendir(dirname)) != NULL)
    {
        rc = 0;
    }
    else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        sfree(gmxdir);
        *p_gmxdir = NULL;
        rc        = EINVAL;
    }
<<<<<<< HEAD
#elif (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64)
    
    if(dirname!=NULL && strlen(dirname)>0)
=======
#elif (defined GMX_NATIVE_WINDOWS)

    if (dirname != NULL && strlen(dirname) > 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        char *     tmpname;
        size_t     namelength;
        int        len;
<<<<<<< HEAD
        
        len = strlen(dirname);
        snew(tmpname,len+3);
        
        strncpy(tmpname,dirname,len+1);
        
        /* Remove possible trailing directory separator */
        if(tmpname[len]=='/' || tmpname[len]=='\\')
        {
            tmpname[len]='\0';
        }
        
        /* Add wildcard */
        strcat(tmpname,"/*");
        
        gmxdir->first = 1;
        if( (gmxdir->windows_handle=_findfirst(tmpname,&gmxdir->finddata))>0L)
=======

        len = strlen(dirname);
        snew(tmpname, len+3);

        strncpy(tmpname, dirname, len+1);

        /* Remove possible trailing directory separator */
        if (tmpname[len] == '/' || tmpname[len] == '\\')
        {
            tmpname[len] = '\0';
        }

        /* Add wildcard */
        strcat(tmpname, "/*");

        gmxdir->first = 1;
        if ( (gmxdir->windows_handle = _findfirst(tmpname, &gmxdir->finddata)) > 0L)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            rc = 0;
        }
        else
        {
<<<<<<< HEAD
            if(errno==EINVAL)
            {
                sfree(gmxdir);
                *p_gmxdir = NULL;
                rc        = EINVAL;                
=======
            if (errno == EINVAL)
            {
                sfree(gmxdir);
                *p_gmxdir = NULL;
                rc        = EINVAL;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            else
            {
                rc        = 0;
            }
        }
    }
    else
    {
        rc = EINVAL;
    }
#else
    gmx_fatal(FARGS,
              "Source compiled without POSIX dirent or windows support - cannot scan directories.\n"
              "In the very unlikely event this is not a compile-time mistake you could consider\n"
              "implementing support for your platform in futil.c, but contact the developers\n"
              "to make sure it's really necessary!\n");
    rc = -1;
#endif
    return rc;
}


int
<<<<<<< HEAD
gmx_directory_nextfile(gmx_directory_t gmxdir,char *name,int maxlength_name)
{
    int                     rc;
    
#ifdef HAVE_DIRENT_H
    
    struct dirent *         direntp_large;
    struct dirent *         p;
    
    
    if(gmxdir!=NULL && gmxdir->dirent_handle!=NULL)
=======
gmx_directory_nextfile(gmx_directory_t gmxdir, char *name, int maxlength_name)
{
    int                     rc;

#ifdef HAVE_DIRENT_H

    struct dirent *         direntp_large;
    struct dirent *         p;


    if (gmxdir != NULL && gmxdir->dirent_handle != NULL)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        /* On some platforms no space is present for d_name in dirent.
         * Since d_name is guaranteed to be the last entry, allocating
         * extra space for dirent will allow more size for d_name.
         * GMX_MAX_PATH should always be >= the max possible d_name.
         */
        smalloc(direntp_large, sizeof(*direntp_large) + GMX_PATH_MAX);
<<<<<<< HEAD
        rc = readdir_r(gmxdir->dirent_handle,direntp_large,&p);

        if(p!=NULL && rc==0)
        {
            strncpy(name,direntp_large->d_name,maxlength_name);
=======
        rc = readdir_r(gmxdir->dirent_handle, direntp_large, &p);

        if (p != NULL && rc == 0)
        {
            strncpy(name, direntp_large->d_name, maxlength_name);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        else
        {
            name[0] = '\0';
            rc      = ENOENT;
        }
        sfree(direntp_large);
    }
<<<<<<< HEAD
    else 
=======
    else
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        name[0] = '\0';
        rc      = EINVAL;
    }
<<<<<<< HEAD
    
#elif (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64)
    
    if(gmxdir!=NULL)
    {
        if(gmxdir->windows_handle<=0)
        {
            
            name[0] = '\0';
            rc      = ENOENT;
        }
        else if(gmxdir->first==1)
        {
            strncpy(name,gmxdir->finddata.name,maxlength_name);
=======

#elif (defined GMX_NATIVE_WINDOWS)

    if (gmxdir != NULL)
    {
        if (gmxdir->windows_handle <= 0)
        {

            name[0] = '\0';
            rc      = ENOENT;
        }
        else if (gmxdir->first == 1)
        {
            strncpy(name, gmxdir->finddata.name, maxlength_name);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            rc            = 0;
            gmxdir->first = 0;
        }
        else
        {
<<<<<<< HEAD
            if(_findnext(gmxdir->windows_handle,&gmxdir->finddata)==0)
            {
                strncpy(name,gmxdir->finddata.name,maxlength_name);
=======
            if (_findnext(gmxdir->windows_handle, &gmxdir->finddata) == 0)
            {
                strncpy(name, gmxdir->finddata.name, maxlength_name);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                rc      = 0;
            }
            else
            {
                name[0] = '\0';
                rc      = ENOENT;
            }
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#else
    gmx_fatal(FARGS,
              "Source compiled without POSIX dirent or windows support - cannot scan directories.\n");
    rc = -1;
#endif
    return rc;
}


<<<<<<< HEAD
int 
=======
int
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
gmx_directory_close(gmx_directory_t gmxdir)
{
    int                     rc;
#ifdef HAVE_DIRENT_H
    rc = (gmxdir != NULL) ? closedir(gmxdir->dirent_handle) : EINVAL;
<<<<<<< HEAD
#elif (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64)
=======
#elif (defined GMX_NATIVE_WINDOWS)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    rc = (gmxdir != NULL) ? _findclose(gmxdir->windows_handle) : EINVAL;
#else
    gmx_fatal(FARGS,
              "Source compiled without POSIX dirent or windows support - cannot scan directories.\n");
    rc = -1;
#endif
<<<<<<< HEAD
    
=======

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    sfree(gmxdir);
    return rc;
}




gmx_bool search_subdirs(const char *parent, char *libdir)
{
<<<<<<< HEAD
    char *ptr;
    gmx_bool found;

    /* Search a few common subdirectory names for the gromacs library dir */
    sprintf(libdir,"%s%cshare%ctop%cgurgle.dat",parent,
            DIR_SEPARATOR,DIR_SEPARATOR,DIR_SEPARATOR);
    found=gmx_fexist(libdir);
    if(!found) {
        sprintf(libdir,"%s%cshare%cgromacs%ctop%cgurgle.dat",parent,
                DIR_SEPARATOR,DIR_SEPARATOR,
                DIR_SEPARATOR,DIR_SEPARATOR);
        found=gmx_fexist(libdir);
    }    
    if(!found) {
        sprintf(libdir,"%s%cshare%cgromacs-%s%ctop%cgurgle.dat",parent,
                DIR_SEPARATOR,DIR_SEPARATOR,VERSION,
                DIR_SEPARATOR,DIR_SEPARATOR);
        found=gmx_fexist(libdir);
    }    
    if(!found) {
        sprintf(libdir,"%s%cshare%cgromacs%cgromacs-%s%ctop%cgurgle.dat",parent,
                DIR_SEPARATOR,DIR_SEPARATOR,DIR_SEPARATOR,
                VERSION,DIR_SEPARATOR,DIR_SEPARATOR);
        found=gmx_fexist(libdir);
    }    

    /* Remove the gurgle.dat part from libdir if we found something */
    if(found) {
        ptr=strrchr(libdir,DIR_SEPARATOR); /* slash or backslash always present, no check necessary */
        *ptr='\0';
=======
    char    *ptr;
    gmx_bool found;

    /* Search a few common subdirectory names for the gromacs library dir */
    sprintf(libdir, "%s%cshare%ctop%cgurgle.dat", parent,
            DIR_SEPARATOR, DIR_SEPARATOR, DIR_SEPARATOR);
    found = gmx_fexist(libdir);
    if (!found)
    {
        sprintf(libdir, "%s%cshare%cgromacs%ctop%cgurgle.dat", parent,
                DIR_SEPARATOR, DIR_SEPARATOR,
                DIR_SEPARATOR, DIR_SEPARATOR);
        found = gmx_fexist(libdir);
    }
    if (!found)
    {
        sprintf(libdir, "%s%cshare%cgromacs-%s%ctop%cgurgle.dat", parent,
                DIR_SEPARATOR, DIR_SEPARATOR, VERSION,
                DIR_SEPARATOR, DIR_SEPARATOR);
        found = gmx_fexist(libdir);
    }
    if (!found)
    {
        sprintf(libdir, "%s%cshare%cgromacs%cgromacs-%s%ctop%cgurgle.dat", parent,
                DIR_SEPARATOR, DIR_SEPARATOR, DIR_SEPARATOR,
                VERSION, DIR_SEPARATOR, DIR_SEPARATOR);
        found = gmx_fexist(libdir);
    }

    /* Remove the gurgle.dat part from libdir if we found something */
    if (found)
    {
        ptr  = strrchr(libdir, DIR_SEPARATOR); /* slash or backslash always present, no check necessary */
        *ptr = '\0';
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    return found;
}


/* Check if the program name begins with "/" on unix/cygwin, or
 * with "\" or "X:\" on windows. If not, the program name
 * is relative to the current directory.
 */
static gmx_bool filename_is_absolute(char *name)
{
<<<<<<< HEAD
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    return ((name[0] == DIR_SEPARATOR) || ((strlen(name)>3) && strncmp(name+1,":\\",2)) == 0);
=======
#ifdef GMX_NATIVE_WINDOWS
    return ((name[0] == DIR_SEPARATOR) || ((strlen(name) > 3) && strncmp(name+1, ":\\", 2)) == 0);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#else
    return (name[0] == DIR_SEPARATOR);
#endif
}

gmx_bool get_libdir(char *libdir)
{
#define GMX_BINNAME_MAX 512
<<<<<<< HEAD
    char bin_name[GMX_BINNAME_MAX];
    char buf[GMX_BINNAME_MAX];
    char full_path[GMX_PATH_MAX+GMX_BINNAME_MAX];
    char system_path[GMX_PATH_MAX];
    char *dir,*ptr,*s,*pdum;
    gmx_bool found=FALSE;
    int i;
=======
    char     bin_name[GMX_BINNAME_MAX];
    char     buf[GMX_BINNAME_MAX];
    char     full_path[GMX_PATH_MAX+GMX_BINNAME_MAX];
    char     system_path[GMX_PATH_MAX];
    char    *dir, *ptr, *s, *pdum;
    gmx_bool found = FALSE;
    int      i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (Program() != NULL)
    {

<<<<<<< HEAD
    /* First - detect binary name */
    if (strlen(Program()) >= GMX_BINNAME_MAX)
    {
        gmx_fatal(FARGS,"The name of the binary is longer than the allowed buffer size (%d):\n'%s'",GMX_BINNAME_MAX,Program());
    }
    strncpy(bin_name,Program(),GMX_BINNAME_MAX-1);

    /* On windows & cygwin we need to add the .exe extension
     * too, or we wont be able to detect that the file exists
     */
#if (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64 || defined __CYGWIN__ || defined __CYGWIN32__)
    if(strlen(bin_name)<3 || gmx_strncasecmp(bin_name+strlen(bin_name)-4,".exe",4))
        strcat(bin_name,".exe");
#endif

    /* Only do the smart search part if we got a real name */
    if (NULL!=bin_name && strncmp(bin_name,"GROMACS",GMX_BINNAME_MAX)) {

        if (!strchr(bin_name,DIR_SEPARATOR)) {
            /* No slash or backslash in name means it must be in the path - search it! */
            /* Add the local dir since it is not in the path on windows */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
            pdum=_getcwd(system_path,sizeof(system_path)-1);
#else
            pdum=getcwd(system_path,sizeof(system_path)-1);
#endif
            sprintf(full_path,"%s%c%s",system_path,DIR_SEPARATOR,bin_name);
            found = gmx_fexist(full_path);
            if (!found && (s=getenv("PATH")) != NULL)
            {
                char *dupped;
                
                dupped=gmx_strdup(s);
                s=dupped;
                while(!found && (dir=gmx_strsep(&s, PATH_SEPARATOR)) != NULL)
                {
                    sprintf(full_path,"%s%c%s",dir,DIR_SEPARATOR,bin_name);
                    found = gmx_fexist(full_path);
                }
                sfree(dupped);
            }
            if (!found)
            {
                return FALSE;
            }
        } else if (!filename_is_absolute(bin_name)) {
            /* name contains directory separators, but 
             * it does not start at the root, i.e.
             * name is relative to the current dir 
             */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
            pdum=_getcwd(buf,sizeof(buf)-1);
#else
            pdum=getcwd(buf,sizeof(buf)-1);
#endif
            sprintf(full_path,"%s%c%s",buf,DIR_SEPARATOR,bin_name);
        } else {
            strncpy(full_path,bin_name,GMX_PATH_MAX);
        }

        /* Now we should have a full path and name in full_path,
         * but on unix it might be a link, or a link to a link to a link..
         */
#if (!defined WIN32 && !defined _WIN32 && !defined WIN64 && !defined _WIN64)
        while( (i=readlink(full_path,buf,sizeof(buf)-1)) > 0 ) {
            buf[i]='\0';
            /* If it doesn't start with "/" it is relative */
            if (buf[0]!=DIR_SEPARATOR) {
                strncpy(strrchr(full_path,DIR_SEPARATOR)+1,buf,GMX_PATH_MAX);
            } else
                strncpy(full_path,buf,GMX_PATH_MAX);
        }
#endif

        /* Remove the executable name - it always contains at least one slash */
        *(strrchr(full_path,DIR_SEPARATOR)+1)='\0';
        /* Now we have the full path to the gromacs executable.
         * Use it to find the library dir. 
         */
        found=FALSE;
        while(!found && ( (ptr=strrchr(full_path,DIR_SEPARATOR)) != NULL ) ) {
            *ptr='\0';
            found=search_subdirs(full_path,libdir);
        }
    }
    }
    /* End of smart searching. If we didn't find it in our parent tree,
     * or if the program name wasn't set, at least try some standard 
     * locations before giving up, in case we are running from e.g. 
     * a users home directory. This only works on unix or cygwin...
     */
#if ((!defined WIN32 && !defined _WIN32 && !defined WIN64 && !defined _WIN64) || defined __CYGWIN__ || defined __CYGWIN32__)
    if(!found) 
        found=search_subdirs("/usr/local",libdir);
    if(!found) 
        found=search_subdirs("/usr",libdir);
    if(!found) 
        found=search_subdirs("/opt",libdir);
=======
        /* First - detect binary name */
        if (strlen(Program()) >= GMX_BINNAME_MAX)
        {
            gmx_fatal(FARGS, "The name of the binary is longer than the allowed buffer size (%d):\n'%s'", GMX_BINNAME_MAX, Program());
        }
        strncpy(bin_name, Program(), GMX_BINNAME_MAX-1);

        /* On windows & cygwin we need to add the .exe extension
         * too, or we wont be able to detect that the file exists
         */
#if (defined GMX_NATIVE_WINDOWS || defined GMX_CYGWIN)
        if (strlen(bin_name) < 3 || gmx_strncasecmp(bin_name+strlen(bin_name)-4, ".exe", 4))
        {
            strcat(bin_name, ".exe");
        }
#endif

        /* Only do the smart search part if we got a real name */
        if (NULL != bin_name && strncmp(bin_name, "GROMACS", GMX_BINNAME_MAX))
        {

            if (!strchr(bin_name, DIR_SEPARATOR))
            {
                /* No slash or backslash in name means it must be in the path - search it! */
                /* Add the local dir since it is not in the path on windows */
#ifdef GMX_NATIVE_WINDOWS
                pdum = _getcwd(system_path, sizeof(system_path)-1);
#else
                pdum = getcwd(system_path, sizeof(system_path)-1);
#endif
                sprintf(full_path, "%s%c%s", system_path, DIR_SEPARATOR, bin_name);
                found = gmx_is_file(full_path);
                if (!found && (s = getenv("PATH")) != NULL)
                {
                    char *dupped;

                    dupped = gmx_strdup(s);
                    s      = dupped;
                    while (!found && (dir = gmx_strsep(&s, PATH_SEPARATOR)) != NULL)
                    {
                        sprintf(full_path, "%s%c%s", dir, DIR_SEPARATOR, bin_name);
                        found = gmx_is_file(full_path);
                    }
                    sfree(dupped);
                }
                if (!found)
                {
                    return FALSE;
                }
            }
            else if (!filename_is_absolute(bin_name))
            {
                /* name contains directory separators, but
                 * it does not start at the root, i.e.
                 * name is relative to the current dir
                 */
#ifdef GMX_NATIVE_WINDOWS
                pdum = _getcwd(buf, sizeof(buf)-1);
#else
                pdum = getcwd(buf, sizeof(buf)-1);
#endif
                sprintf(full_path, "%s%c%s", buf, DIR_SEPARATOR, bin_name);
            }
            else
            {
                strncpy(full_path, bin_name, GMX_PATH_MAX);
            }

            /* Now we should have a full path and name in full_path,
             * but on unix it might be a link, or a link to a link to a link..
             */
#ifndef GMX_NATIVE_WINDOWS
            while ( (i = readlink(full_path, buf, sizeof(buf)-1)) > 0)
            {
                buf[i] = '\0';
                /* If it doesn't start with "/" it is relative */
                if (buf[0] != DIR_SEPARATOR)
                {
                    strncpy(strrchr(full_path, DIR_SEPARATOR)+1, buf, GMX_PATH_MAX);
                }
                else
                {
                    strncpy(full_path, buf, GMX_PATH_MAX);
                }
            }
#endif

            /* Remove the executable name - it always contains at least one slash */
            *(strrchr(full_path, DIR_SEPARATOR)+1) = '\0';
            /* Now we have the full path to the gromacs executable.
             * Use it to find the library dir.
             */
            found = FALSE;
            while (!found && ( (ptr = strrchr(full_path, DIR_SEPARATOR)) != NULL ) )
            {
                *ptr  = '\0';
                found = search_subdirs(full_path, libdir);
            }
        }
    }
    /* End of smart searching. If we didn't find it in our parent tree,
     * or if the program name wasn't set, at least try some standard
     * locations before giving up, in case we are running from e.g.
     * a users home directory. This only works on unix or cygwin...
     */
#ifndef GMX_NATIVE_WINDOWS
    if (!found)
    {
        found = search_subdirs("/usr/local", libdir);
    }
    if (!found)
    {
        found = search_subdirs("/usr", libdir);
    }
    if (!found)
    {
        found = search_subdirs("/opt", libdir);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
    return found;
}


char *low_gmxlibfn(const char *file, gmx_bool bAddCWD, gmx_bool bFatal)
{
<<<<<<< HEAD
    char *ret;
    char *lib,*dir;
    char buf[1024];
    char libpath[GMX_PATH_MAX];
    gmx_bool env_is_set=FALSE;
    char   *s,tmppath[GMX_PATH_MAX];

    /* GMXLIB can be a path now */
    lib=getenv("GMXLIB");
    if (lib != NULL)
    {
        env_is_set=TRUE;
        strncpy(libpath,lib,GMX_PATH_MAX);
    } 
    else if (!get_libdir(libpath))
    {
        strncpy(libpath,GMXLIBDIR,GMX_PATH_MAX);
=======
    char    *ret;
    char    *lib, *dir;
    char     buf[1024];
    char     libpath[GMX_PATH_MAX];
    gmx_bool env_is_set = FALSE;
    char    *s, tmppath[GMX_PATH_MAX];

    /* GMXLIB can be a path now */
    lib = getenv("GMXLIB");
    if (lib != NULL)
    {
        env_is_set = TRUE;
        strncpy(libpath, lib, GMX_PATH_MAX);
    }
    else if (!get_libdir(libpath))
    {
        strncpy(libpath, GMXLIBDIR, GMX_PATH_MAX);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    ret = NULL;
    if (bAddCWD && gmx_fexist(file))
    {
        ret = gmx_strdup(file);
    }
<<<<<<< HEAD
    else 
    {
        strncpy(tmppath,libpath,GMX_PATH_MAX);
        s=tmppath;
        while(ret == NULL && (dir=gmx_strsep(&s, PATH_SEPARATOR)) != NULL )
        {
            sprintf(buf,"%s%c%s",dir,DIR_SEPARATOR,file);
=======
    else
    {
        strncpy(tmppath, libpath, GMX_PATH_MAX);
        s = tmppath;
        while (ret == NULL && (dir = gmx_strsep(&s, PATH_SEPARATOR)) != NULL)
        {
            sprintf(buf, "%s%c%s", dir, DIR_SEPARATOR, file);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            if (gmx_fexist(buf))
            {
                ret = gmx_strdup(buf);
            }
        }
<<<<<<< HEAD
        if (ret == NULL && bFatal) 
        {
            if (env_is_set) 
=======
        if (ret == NULL && bFatal)
        {
            if (env_is_set)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                gmx_fatal(FARGS,
                          "Library file %s not found %sin your GMXLIB path.",
                          file, bAddCWD ? "in current dir nor " : "");
            }
            else
            {
                gmx_fatal(FARGS,
                          "Library file %s not found %sin default directories.\n"
<<<<<<< HEAD
                        "(You can set the directories to search with the GMXLIB path variable)",
=======
                          "(You can set the directories to search with the GMXLIB path variable)",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                          file, bAddCWD ? "in current dir nor " : "");
            }
        }
    }

    return ret;
}





<<<<<<< HEAD
FILE *low_libopen(const char *file,gmx_bool bFatal)
=======
FILE *low_libopen(const char *file, gmx_bool bFatal)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
    FILE *ff;
    char *fn;

<<<<<<< HEAD
    fn=low_gmxlibfn(file,TRUE,bFatal);

    if (fn==NULL) {
        ff=NULL;
    } else {
      if (debug)
	fprintf(debug,"Opening library file %s\n",fn);
      ff=fopen(fn,"r");
=======
    fn = low_gmxlibfn(file, TRUE, bFatal);

    if (fn == NULL)
    {
        ff = NULL;
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "Opening library file %s\n", fn);
        }
        ff = fopen(fn, "r");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    sfree(fn);

    return ff;
}

char *gmxlibfn(const char *file)
{
<<<<<<< HEAD
    return low_gmxlibfn(file,TRUE,TRUE);
=======
    return low_gmxlibfn(file, TRUE, TRUE);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

FILE *libopen(const char *file)
{
<<<<<<< HEAD
    return low_libopen(file,TRUE);
=======
    return low_libopen(file, TRUE);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

void gmx_tmpnam(char *buf)
{
<<<<<<< HEAD
    int i,len,fd;

    if ((len = strlen(buf)) < 7)
        gmx_fatal(FARGS,"Buf passed to gmx_tmpnam must be at least 7 bytes long");
    for(i=len-6; (i<len); i++) {
        buf[i] = 'X';
    }
    /* mktemp is dangerous and we should use mkstemp instead, but 
     * since windows doesnt support it we have to separate the cases.
     * 20090307: mktemp deprecated, use iso c++ _mktemp instead.
     */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
=======
    int i, len, fd;

    if ((len = strlen(buf)) < 7)
    {
        gmx_fatal(FARGS, "Buf passed to gmx_tmpnam must be at least 7 bytes long");
    }
    for (i = len-6; (i < len); i++)
    {
        buf[i] = 'X';
    }
    /* mktemp is dangerous and we should use mkstemp instead, but
     * since windows doesnt support it we have to separate the cases.
     * 20090307: mktemp deprecated, use iso c++ _mktemp instead.
     */
#ifdef GMX_NATIVE_WINDOWS
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    _mktemp(buf);
#else
    fd = mkstemp(buf);

<<<<<<< HEAD
    switch (fd) {
        case EINVAL:
            gmx_fatal(FARGS,"Invalid template %s for mkstemp",buf);
            break;
        case EEXIST:
            gmx_fatal(FARGS,"mkstemp created existing file",buf);
            break;
        case EACCES: 
            gmx_fatal(FARGS,"Permission denied for opening %s",buf);
            break;
        default:
            break;
    }   
=======
    switch (fd)
    {
        case EINVAL:
            gmx_fatal(FARGS, "Invalid template %s for mkstemp", buf);
            break;
        case EEXIST:
            gmx_fatal(FARGS, "mkstemp created existing file", buf);
            break;
        case EACCES:
            gmx_fatal(FARGS, "Permission denied for opening %s", buf);
            break;
        default:
            break;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    close(fd);
#endif
    /* name in Buf should now be OK */
}

int gmx_truncatefile(char *path, gmx_off_t length)
{
#ifdef _MSC_VER
    /* Microsoft visual studio does not have "truncate" */
<<<<<<< HEAD
    HANDLE fh;
=======
    HANDLE        fh;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    LARGE_INTEGER win_length;

    win_length.QuadPart = length;

<<<<<<< HEAD
    fh = CreateFile(path,GENERIC_READ | GENERIC_WRITE,0,NULL,
            OPEN_EXISTING,0,NULL);
    SetFilePointerEx(fh,win_length,NULL,FILE_BEGIN);
=======
    fh = CreateFile(path, GENERIC_READ | GENERIC_WRITE, 0, NULL,
                    OPEN_EXISTING, 0, NULL);
    SetFilePointerEx(fh, win_length, NULL, FILE_BEGIN);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    SetEndOfFile(fh);
    CloseHandle(fh);

    return 0;
#else
<<<<<<< HEAD
    return truncate(path,length);
=======
    return truncate(path, length);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}


int gmx_file_rename(const char *oldname, const char *newname)
{
<<<<<<< HEAD
#if (!(defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)))
    /* under unix, rename() is atomic (at least, it should be). */
    return rename(oldname, newname);
#else
    if (MoveFileEx(oldname, newname, 
                   MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
        return 0;
    else
        return 1;
=======
#ifndef GMX_NATIVE_WINDOWS
    /* under unix, rename() is atomic (at least, it should be). */
    return rename(oldname, newname);
#else
    if (MoveFileEx(oldname, newname,
                   MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
    {
        return 0;
    }
    else
    {
        return 1;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

int gmx_file_copy(const char *oldname, const char *newname, gmx_bool copy_if_empty)
{
/* the full copy buffer size: */
#define FILECOPY_BUFSIZE (1<<16)
<<<<<<< HEAD
    FILE *in=NULL; 
    FILE *out=NULL;
    char *buf;

    snew(buf, FILECOPY_BUFSIZE); 

    in=fopen(oldname, "rb");
    if (!in)
        goto error;
=======
    FILE *in  = NULL;
    FILE *out = NULL;
    char *buf;

    snew(buf, FILECOPY_BUFSIZE);

    in = fopen(oldname, "rb");
    if (!in)
    {
        goto error;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* If we don't copy when empty, we postpone opening the file
       until we're actually ready to write. */
    if (copy_if_empty)
    {
<<<<<<< HEAD
        out=fopen(newname, "wb");
        if (!out)
            goto error;
    }

    while(!feof(in))
    {
        size_t nread;
        
        nread=fread(buf, sizeof(char), FILECOPY_BUFSIZE, in);
        if (nread>0)
=======
        out = fopen(newname, "wb");
        if (!out)
        {
            goto error;
        }
    }

    while (!feof(in))
    {
        size_t nread;

        nread = fread(buf, sizeof(char), FILECOPY_BUFSIZE, in);
        if (nread > 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            size_t ret;
            if (!out)
            {
                /* so this is where we open when copy_if_empty is false:
                   here we know we read something. */
<<<<<<< HEAD
                out=fopen(newname, "wb");
                if (!out)
                    goto error;
            }
            ret=fwrite(buf, sizeof(char), nread, out);
            if (ret!=nread)
=======
                out = fopen(newname, "wb");
                if (!out)
                {
                    goto error;
                }
            }
            ret = fwrite(buf, sizeof(char), nread, out);
            if (ret != nread)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                goto error;
            }
        }
        if (ferror(in))
<<<<<<< HEAD
            goto error;
=======
        {
            goto error;
        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    sfree(buf);
    fclose(in);
    fclose(out);
    return 0;
error:
    sfree(buf);
    if (in)
<<<<<<< HEAD
        fclose(in);
    if (out)
        fclose(out);
=======
    {
        fclose(in);
    }
    if (out)
    {
        fclose(out);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    return 1;
#undef FILECOPY_BUFSIZE
}


int gmx_fsync(FILE *fp)
{
<<<<<<< HEAD
    int rc=0;

#ifdef GMX_FAHCORE
    /* the fahcore defines its own os-independent fsync */
    rc=fah_fsync(fp);
#else /* GMX_FAHCORE */
    {
        int fn=-1;

        /* get the file number */
#if defined(HAVE_FILENO)
        fn= fileno(fp);
#elif defined(HAVE__FILENO)
        fn= _fileno(fp);
=======
    int rc = 0;

#ifdef GMX_FAHCORE
    /* the fahcore defines its own os-independent fsync */
    rc = fah_fsync(fp);
#else /* GMX_FAHCORE */
    {
        int fn = -1;

        /* get the file number */
#if defined(HAVE_FILENO)
        fn = fileno(fp);
#elif defined(HAVE__FILENO)
        fn = _fileno(fp);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

        /* do the actual fsync */
        if (fn >= 0)
        {
#if (defined(HAVE_FSYNC))
<<<<<<< HEAD
            rc=fsync(fn);
#elif (defined(HAVE__COMMIT)) 
            rc=_commit(fn);
=======
            rc = fsync(fn);
#elif (defined(HAVE__COMMIT))
            rc = _commit(fn);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
        }
    }
#endif /* GMX_FAHCORE */

    /* We check for these error codes this way because POSIX requires them
       to be defined, and using anything other than macros is unlikely: */
#ifdef EINTR
    /* we don't want to report an error just because fsync() caught a signal.
       For our purposes, we can just ignore this. */
<<<<<<< HEAD
    if (rc && errno==EINTR)
        rc=0;
#endif
#ifdef EINVAL
    /* we don't want to report an error just because we tried to fsync() 
       stdout, a socket or a pipe. */
    if (rc && errno==EINVAL)
        rc=0;
#endif
    return rc;
}



=======
    if (rc && errno == EINTR)
    {
        rc = 0;
    }
#endif
#ifdef EINVAL
    /* we don't want to report an error just because we tried to fsync()
       stdout, a socket or a pipe. */
    if (rc && errno == EINVAL)
    {
        rc = 0;
    }
#endif
    return rc;
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
