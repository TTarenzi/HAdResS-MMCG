<<<<<<< HEAD
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
 */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"

void 
F77_FUNC(slasd6,SLASD6)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	float *d__, 
	float *vf, 
	float *vl, 
	float *alpha, 
	float *beta, 
	int *idxq, 
	int *perm, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	float *givnum,
	int *ldgnum, 
	float *poles, 
	float *difl, 
	float *difr, 
	float *z__, 
	int *k, 
	float *c__, 
	float *s, 
	float *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, i__1;
    float d__1, d__2;

    int i__, m, n, n1, n2, iw, idx, idxc, idxp, ivfw, ivlw;
    int isigma;
    float orgnrm;
    int c__0 = 0;
    float one = 1.0;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    --vf;
    --vl;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --difr;
    --z__;
    --work;
    --iwork;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;

    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;

    d__1 = fabs(*alpha); 
    d__2 = fabs(*beta);
    orgnrm = (d__1 > d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      d__1 = fabs(d__[i__]);
	if (d__1 > orgnrm)
	    orgnrm = d__1;
    }
    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    F77_FUNC(slasd7,SLASD7)(icompq, nl, nr, sqre, k, &d__[1], &z__[1], &work[iw], &vf[1], &
	    work[ivfw], &vl[1], &work[ivlw], alpha, beta, &work[isigma], &
	    iwork[idx], &iwork[idxp], &idxq[1], &perm[1], givptr, &givcol[
	    givcol_offset], ldgcol, &givnum[givnum_offset], ldgnum, c__, s, 
	    info);

    F77_FUNC(slasd8,SLASD8)(icompq, k, &d__[1], &z__[1], &vf[1], &vl[1], &difl[1], &difr[1], 
	    ldgnum, &work[isigma], &work[iw], info);

    if (*icompq == 1) {
	F77_FUNC(scopy,SCOPY)(k, &d__[1], &c__1, &poles[poles_dim1 + 1], &c__1);
	F77_FUNC(scopy,SCOPY)(k, &work[isigma], &c__1, &poles[(poles_dim1 << 1) + 1], &c__1);
    }

    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = *k;
    n2 = n - *k;
    F77_FUNC(slamrg,SLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}


