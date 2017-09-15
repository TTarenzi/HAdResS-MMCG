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
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
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

#ifndef _domdec_network_h
#define _domdec_network_h
<<<<<<< HEAD

#include "typedefs.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif
=======
#include "visibility.h"
#include "typedefs.h"
#include "types/commrec.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
extern "C" {
#endif

enum {
<<<<<<< HEAD
    dddirForward,dddirBackward
=======
    dddirForward, dddirBackward
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
};

/* Move integers in the comm. region one cell along the domain decomposition
 * in the dimension indexed by ddimind
 * forward (direction=dddirFoward) or backward (direction=dddirBackward).
 */
void
dd_sendrecv_int(const gmx_domdec_t *dd,
<<<<<<< HEAD
                int ddimind,int direction,
                int *buf_s,int n_s,
                int *buf_r,int n_r);
=======
                int ddimind, int direction,
                int *buf_s, int n_s,
                int *buf_r, int n_r);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/* Move reals in the comm. region one cell along the domain decomposition
 * in the dimension indexed by ddimind
 * forward (direction=dddirFoward) or backward (direction=dddirBackward).
 */
void
dd_sendrecv_real(const gmx_domdec_t *dd,
<<<<<<< HEAD
                 int ddimind,int direction,
                 real *buf_s,int n_s,
                 real *buf_r,int n_r);
=======
                 int ddimind, int direction,
                 real *buf_s, int n_s,
                 real *buf_r, int n_r);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/* Move revc's in the comm. region one cell along the domain decomposition
 * in dimension indexed by ddimind
 * forward (direction=dddirFoward) or backward (direction=dddirBackward).
 */
void
dd_sendrecv_rvec(const gmx_domdec_t *dd,
<<<<<<< HEAD
                 int ddimind,int direction,
                 rvec *buf_s,int n_s,
                 rvec *buf_r,int n_r);
=======
                 int ddimind, int direction,
                 rvec *buf_s, int n_s,
                 rvec *buf_r, int n_r);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/* Move revc's in the comm. region one cell along the domain decomposition
 * in dimension indexed by ddimind
 * simultaneously in the forward and backward directions.
 */
void
dd_sendrecv2_rvec(const gmx_domdec_t *dd,
<<<<<<< HEAD
		  int ddimind,
		  rvec *buf_s_fw,int n_s_fw,
		  rvec *buf_r_fw,int n_r_fw,
		  rvec *buf_s_bw,int n_s_bw,
		  rvec *buf_r_bw,int n_r_bw);
=======
                  int ddimind,
                  rvec *buf_s_fw, int n_s_fw,
                  rvec *buf_r_fw, int n_r_fw,
                  rvec *buf_s_bw, int n_s_bw,
                  rvec *buf_r_bw, int n_r_bw);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/* The functions below perform the same operations as the MPI functions
 * with the same name appendices, but over the domain decomposition
 * nodes only.
 * The DD master node is the master for these operations.
 */

<<<<<<< HEAD
void
dd_bcast(gmx_domdec_t *dd,int nbytes,void *data);

/* Copies src to dest on the master node and then broadcasts */
void
dd_bcastc(gmx_domdec_t *dd,int nbytes,void *src,void *dest);

void
dd_scatter(gmx_domdec_t *dd,int nbytes,void *src,void *dest);

void
dd_gather(gmx_domdec_t *dd,int nbytes,void *src,void *dest);
=======
GMX_LIBMD_EXPORT
void
dd_bcast(gmx_domdec_t *dd, int nbytes, void *data);

/* Copies src to dest on the master node and then broadcasts */
void
dd_bcastc(gmx_domdec_t *dd, int nbytes, void *src, void *dest);

void
dd_scatter(gmx_domdec_t *dd, int nbytes, void *src, void *dest);

void
dd_gather(gmx_domdec_t *dd, int nbytes, void *src, void *dest);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/* If rcount==0, rbuf is allowed to be NULL */
void
dd_scatterv(gmx_domdec_t *dd,
<<<<<<< HEAD
            int *scounts,int *disps,void *sbuf,
            int rcount,void *rbuf);
=======
            int *scounts, int *disps, void *sbuf,
            int rcount, void *rbuf);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

/* If scount==0, sbuf is allowed to be NULL */
void
dd_gatherv(gmx_domdec_t *dd,
<<<<<<< HEAD
           int scount,void *sbuf,
           int *rcounts,int *disps,void *rbuf);
=======
           int scount, void *sbuf,
           int *rcounts, int *disps, void *rbuf);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _domdec_network_h */
=======
#endif  /* _domdec_network_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
