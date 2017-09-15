<<<<<<< HEAD
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2008
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

#include <string.h>
#include "domdec_network.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "tmpi.h"
#endif


#define DDMASTERRANK(dd)   (dd->masterrank)


void dd_sendrecv_int(const gmx_domdec_t *dd,
<<<<<<< HEAD
                     int ddimind,int direction,
                     int *buf_s,int n_s,
                     int *buf_r,int n_r)
{
#ifdef GMX_MPI
    int rank_s,rank_r;
    MPI_Status stat;
    
    rank_s = dd->neighbor[ddimind][direction==dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction==dddirForward ? 1 : 0];
    
    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s,n_s*sizeof(int),MPI_BYTE,rank_s,0,
                     buf_r,n_r*sizeof(int),MPI_BYTE,rank_r,0,
                     dd->mpi_comm_all,&stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s,n_s*sizeof(int),MPI_BYTE,rank_s,0,
=======
                     int ddimind, int direction,
                     int *buf_s, int n_s,
                     int *buf_r, int n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
<<<<<<< HEAD
        MPI_Recv(    buf_r,n_r*sizeof(int),MPI_BYTE,rank_r,0,
                     dd->mpi_comm_all,&stat);
    }
    
=======
        MPI_Recv(    buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

void dd_sendrecv_real(const gmx_domdec_t *dd,
<<<<<<< HEAD
                      int ddimind,int direction,
                      real *buf_s,int n_s,
                      real *buf_r,int n_r)
{
#ifdef GMX_MPI
    int rank_s,rank_r;
    MPI_Status stat;
    
    rank_s = dd->neighbor[ddimind][direction==dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction==dddirForward ? 1 : 0];
    
    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s,n_s*sizeof(real),MPI_BYTE,rank_s,0,
                     buf_r,n_r*sizeof(real),MPI_BYTE,rank_r,0,
                     dd->mpi_comm_all,&stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s,n_s*sizeof(real),MPI_BYTE,rank_s,0,
=======
                      int ddimind, int direction,
                      real *buf_s, int n_s,
                      real *buf_r, int n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
<<<<<<< HEAD
        MPI_Recv(    buf_r,n_r*sizeof(real),MPI_BYTE,rank_r,0,
                     dd->mpi_comm_all,&stat);
    }
    
=======
        MPI_Recv(    buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

void dd_sendrecv_rvec(const gmx_domdec_t *dd,
<<<<<<< HEAD
                      int ddimind,int direction,
                      rvec *buf_s,int n_s,
                      rvec *buf_r,int n_r)
{
#ifdef GMX_MPI
    int rank_s,rank_r;
    MPI_Status stat;
    
    rank_s = dd->neighbor[ddimind][direction==dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction==dddirForward ? 1 : 0];
    
    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s[0],n_s*sizeof(rvec),MPI_BYTE,rank_s,0,
                     buf_r[0],n_r*sizeof(rvec),MPI_BYTE,rank_r,0,
                     dd->mpi_comm_all,&stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s[0],n_s*sizeof(rvec),MPI_BYTE,rank_s,0,
                     dd->mpi_comm_all);
    } else if (n_r)
    {
        MPI_Recv(    buf_r[0],n_r*sizeof(rvec),MPI_BYTE,rank_r,0,
                     dd->mpi_comm_all,&stat);
    }
    
=======
                      int ddimind, int direction,
                      rvec *buf_s, int n_s,
                      rvec *buf_r, int n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s[0], n_s*sizeof(rvec), MPI_BYTE, rank_s, 0,
                     buf_r[0], n_r*sizeof(rvec), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s[0], n_s*sizeof(rvec), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r[0], n_r*sizeof(rvec), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

void dd_sendrecv2_rvec(const gmx_domdec_t *dd,
                       int ddimind,
<<<<<<< HEAD
                       rvec *buf_s_fw,int n_s_fw,
                       rvec *buf_r_fw,int n_r_fw,
                       rvec *buf_s_bw,int n_s_bw,
                       rvec *buf_r_bw,int n_r_bw)
{
#ifdef GMX_MPI
    int rank_fw,rank_bw,nreq;
    MPI_Request req[4];
    MPI_Status  stat[4];
    
    rank_fw = dd->neighbor[ddimind][0];
    rank_bw = dd->neighbor[ddimind][1];
    
=======
                       rvec *buf_s_fw, int n_s_fw,
                       rvec *buf_r_fw, int n_r_fw,
                       rvec *buf_s_bw, int n_s_bw,
                       rvec *buf_r_bw, int n_r_bw)
{
#ifdef GMX_MPI
    int         rank_fw, rank_bw, nreq;
    MPI_Request req[4];
    MPI_Status  stat[4];

    rank_fw = dd->neighbor[ddimind][0];
    rank_bw = dd->neighbor[ddimind][1];

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if (!dd->bSendRecv2)
    {
        /* Try to send and receive in two directions simultaneously.
         * Should be faster, especially on machines
         * with full 3D communication networks.
         * However, it could be that communication libraries are
         * optimized for MPI_Sendrecv and non-blocking MPI calls
         * are slower.
         * SendRecv2 can be turned on with the env.var. GMX_DD_SENDRECV2
         */
        nreq = 0;
        if (n_r_fw)
        {
<<<<<<< HEAD
            MPI_Irecv(buf_r_fw[0],n_r_fw*sizeof(rvec),MPI_BYTE,
                      rank_bw,0,dd->mpi_comm_all,&req[nreq++]);
        }
        if (n_r_bw)
        {
            MPI_Irecv(buf_r_bw[0],n_r_bw*sizeof(rvec),MPI_BYTE,
                      rank_fw,1,dd->mpi_comm_all,&req[nreq++]);
        }
        if (n_s_fw)
        {
            MPI_Isend(buf_s_fw[0],n_s_fw*sizeof(rvec),MPI_BYTE,
                      rank_fw,0,dd->mpi_comm_all,&req[nreq++]);
        }
        if (n_s_bw)
        {
            MPI_Isend(buf_s_bw[0],n_s_bw*sizeof(rvec),MPI_BYTE,
                      rank_bw,1,dd->mpi_comm_all,&req[nreq++]);
        }
        if (nreq)
        {
            MPI_Waitall(nreq,req,stat);
=======
            MPI_Irecv(buf_r_fw[0], n_r_fw*sizeof(rvec), MPI_BYTE,
                      rank_bw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_r_bw)
        {
            MPI_Irecv(buf_r_bw[0], n_r_bw*sizeof(rvec), MPI_BYTE,
                      rank_fw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_fw)
        {
            MPI_Isend(buf_s_fw[0], n_s_fw*sizeof(rvec), MPI_BYTE,
                      rank_fw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_bw)
        {
            MPI_Isend(buf_s_bw[0], n_s_bw*sizeof(rvec), MPI_BYTE,
                      rank_bw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (nreq)
        {
            MPI_Waitall(nreq, req, stat);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
    else
    {
        /* Communicate in two ordered phases.
         * This is slower, even on a dual-core Opteron cluster
         * with a single full-duplex network connection per machine.
         */
        /* Forward */
<<<<<<< HEAD
        MPI_Sendrecv(buf_s_fw[0],n_s_fw*sizeof(rvec),MPI_BYTE,rank_fw,0,
                     buf_r_fw[0],n_r_fw*sizeof(rvec),MPI_BYTE,rank_bw,0,
                     dd->mpi_comm_all,&stat[0]);
        /* Backward */
        MPI_Sendrecv(buf_s_bw[0],n_s_bw*sizeof(rvec),MPI_BYTE,rank_bw,0,
                     buf_r_bw[0],n_r_bw*sizeof(rvec),MPI_BYTE,rank_fw,0,
                     dd->mpi_comm_all,&stat[0]);
=======
        MPI_Sendrecv(buf_s_fw[0], n_s_fw*sizeof(rvec), MPI_BYTE, rank_fw, 0,
                     buf_r_fw[0], n_r_fw*sizeof(rvec), MPI_BYTE, rank_bw, 0,
                     dd->mpi_comm_all, &stat[0]);
        /* Backward */
        MPI_Sendrecv(buf_s_bw[0], n_s_bw*sizeof(rvec), MPI_BYTE, rank_bw, 0,
                     buf_r_bw[0], n_r_bw*sizeof(rvec), MPI_BYTE, rank_fw, 0,
                     dd->mpi_comm_all, &stat[0]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
#endif
}

/* IBM's BlueGene(/L) MPI_Bcast dereferences the data pointer
 * even when 0 == nbytes, so we protect calls to it on BlueGene.
 * Fortunately dd_bcast() and dd_bcastc() are only
 * called during DD setup and partition.
 */

<<<<<<< HEAD
void dd_bcast(gmx_domdec_t *dd,int nbytes,void *data)
=======
void dd_bcast(gmx_domdec_t *dd, int nbytes, void *data)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
#ifdef GMX_MPI
#ifdef GMX_BLUEGENE
    if (nbytes > 0)
    {
#endif
<<<<<<< HEAD
    MPI_Bcast(data,nbytes,MPI_BYTE,
              DDMASTERRANK(dd),dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
    }
=======
    MPI_Bcast(data, nbytes, MPI_BYTE,
              DDMASTERRANK(dd), dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
#endif
}

<<<<<<< HEAD
void dd_bcastc(gmx_domdec_t *dd,int nbytes,void *src,void *dest)
{
    if (DDMASTER(dd))
    {
        memcpy(dest,src,nbytes);
=======
void dd_bcastc(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
    if (DDMASTER(dd))
    {
        memcpy(dest, src, nbytes);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
#ifdef GMX_MPI
#ifdef GMX_BLUEGENE
    if (nbytes > 0)
    {
#endif
<<<<<<< HEAD
    MPI_Bcast(dest,nbytes,MPI_BYTE,
              DDMASTERRANK(dd),dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
    }
=======
    MPI_Bcast(dest, nbytes, MPI_BYTE,
              DDMASTERRANK(dd), dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
#endif
}

<<<<<<< HEAD
void dd_scatter(gmx_domdec_t *dd,int nbytes,void *src,void *dest)
{
#ifdef GMX_MPI
    MPI_Scatter(src,nbytes,MPI_BYTE,
                dest,nbytes,MPI_BYTE,
                DDMASTERRANK(dd),dd->mpi_comm_all);
#endif
}

void dd_gather(gmx_domdec_t *dd,int nbytes,void *src,void *dest)
{
#ifdef GMX_MPI
    MPI_Gather(src,nbytes,MPI_BYTE,
               dest,nbytes,MPI_BYTE,
               DDMASTERRANK(dd),dd->mpi_comm_all);
=======
void dd_scatter(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
#ifdef GMX_MPI
    MPI_Scatter(src, nbytes, MPI_BYTE,
                dest, nbytes, MPI_BYTE,
                DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}

void dd_gather(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
#ifdef GMX_MPI
    MPI_Gather(src, nbytes, MPI_BYTE,
               dest, nbytes, MPI_BYTE,
               DDMASTERRANK(dd), dd->mpi_comm_all);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

void dd_scatterv(gmx_domdec_t *dd,
<<<<<<< HEAD
                 int *scounts,int *disps,void *sbuf,
                 int rcount,void *rbuf)
=======
                 int *scounts, int *disps, void *sbuf,
                 int rcount, void *rbuf)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
#ifdef GMX_MPI
    int dum;

    if (rcount == 0)
    {
        /* MPI does not allow NULL pointers */
        rbuf = &dum;
    }
<<<<<<< HEAD
    MPI_Scatterv(sbuf,scounts,disps,MPI_BYTE,
                 rbuf,rcount,MPI_BYTE,
                 DDMASTERRANK(dd),dd->mpi_comm_all);
=======
    MPI_Scatterv(sbuf, scounts, disps, MPI_BYTE,
                 rbuf, rcount, MPI_BYTE,
                 DDMASTERRANK(dd), dd->mpi_comm_all);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

void dd_gatherv(gmx_domdec_t *dd,
<<<<<<< HEAD
                int scount,void *sbuf,
                int *rcounts,int *disps,void *rbuf)
=======
                int scount, void *sbuf,
                int *rcounts, int *disps, void *rbuf)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
{
#ifdef GMX_MPI
    int dum;

    if (scount == 0)
    {
        /* MPI does not allow NULL pointers */
        sbuf = &dum;
    }
<<<<<<< HEAD
    MPI_Gatherv(sbuf,scount,MPI_BYTE,
                rbuf,rcounts,disps,MPI_BYTE,
                DDMASTERRANK(dd),dd->mpi_comm_all);
#endif
}

=======
    MPI_Gatherv(sbuf, scount, MPI_BYTE,
                rbuf, rcounts, disps, MPI_BYTE,
                DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
