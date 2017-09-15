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
 *                        VERSION 4.0.5
 * Written by Christoph Junghans, Brad Lambeth, and possibly others.
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * All rights reserved.

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
 */
 
=======
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
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
#include "adress.h"
#include "maths.h"
#include "pbc.h"
#include "types/simple.h"
#include "typedefs.h"
#include "vec.h"

<<<<<<< HEAD
real 
adress_weight(rvec            x,
              int             adresstype,
              real            adressr,
              real            adressw,
              rvec *          ref,
              t_pbc *         pbc,
=======
real
adress_weight(rvec                 x,
              int                  adresstype,
              real                 adressr,
              real                 adressw,
              rvec      *          ref,
              t_pbc      *         pbc,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
              t_forcerec *         fr )
{
    int  i;
    real l2 = adressr+adressw;
<<<<<<< HEAD
    real sqr_dl,dl;
    real tmp;
    rvec dx;
    real H5, H, H2, l;

    sqr_dl = 0.0;

    if (pbc) 
    {
        pbc_dx(pbc,(*ref),x,dx);
    } 
    else 
    {
        rvec_sub((*ref),x,dx);
    }

    switch(adresstype)
    {
    case eAdressOff:
        /* default to explicit simulation */
        return 1;
    case eAdressConst:              
        /* constant value for weighting function = adressw */
        return fr->adress_const_wf;
    case eAdressXSplit:              
        /* plane through center of ref, varies in x direction */
        sqr_dl         = dx[0]*dx[0];
        break;
    case eAdressSphere:
        /* point at center of ref, assuming cubic geometry */
        for(i=0;i<3;i++){
            sqr_dl    += dx[i]*dx[i];
        }
        break;
    default:
        /* default to explicit simulation */
        return 1;
    }
    
    dl=sqrt(sqr_dl);
    
    /* molecule is coarse grained */
    if (dl > l2)
    {
        return 0;
    }
    /* molecule is explicit */
    else if (dl < adressr)
    {
        return 1;
    }
    /* hybrid region */
    else
    {
        //tmp=cos((dl-adressr)*M_PI/2/adressw);
        //return tmp*tmp;
        l=fabs(dl-adressr);
        H=adressw;
        H5=H*H*H*H*H;
        
        tmp=(1-30.0/H5*(l*l*l*l*l/5.0-1.0/2.0*l*l*l*l*H+l*l*l/3.0*H*H));
        //printf ("DD %g %g\n", l, tmp);
        return tmp;
        
        
    }
}

real
Dadress_weight(rvec            x,
              int             adresstype,
              real            adressr,
              real            adressw,
              rvec *          ref,
              t_pbc *         pbc,
              t_forcerec *         fr )
{
        int  i;
    real l2 = adressr+adressw;
    real sqr_dl,dl;
    real tmp;
    rvec dx;
    real H5, l, H;
    
    sqr_dl = 0.0;

    if (pbc)
    {
        pbc_dx(pbc,(*ref),x,dx);
    }
    else
    {
        rvec_sub((*ref),x,dx);
    }

    switch(adresstype)
    {
    case eAdressOff:
        /* default to explicit simulation */
        return 0;
    case eAdressConst:
        /* constant value for weighting function = adressw */
        return fr->adress_const_wf;
    case eAdressXSplit:
        /* plane through center of ref, varies in x direction */
        sqr_dl         = dx[0]*dx[0];
        break;
    case eAdressSphere:
        /* point at center of ref, assuming cubic geometry */
        for(i=0;i<3;i++){
            sqr_dl    += dx[i]*dx[i];
        }
        break;
    default:
        /* default to explicit simulation */
        return 0;
    }

    dl=sqrt(sqr_dl);
=======
    real sqr_dl, dl;
    real tmp;
    rvec dx;

    sqr_dl = 0.0;

    if (pbc)
    {
        pbc_dx(pbc, (*ref), x, dx);
    }
    else
    {
        rvec_sub((*ref), x, dx);
    }

    switch (adresstype)
    {
        case eAdressOff:
            /* default to explicit simulation */
            return 1;
        case eAdressConst:
            /* constant value for weighting function = adressw */
            return fr->adress_const_wf;
        case eAdressXSplit:
            /* plane through center of ref, varies in x direction */
            sqr_dl         = dx[0]*dx[0];
            break;
        case eAdressSphere:
            /* point at center of ref, assuming cubic geometry */
            for (i = 0; i < 3; i++)
            {
                sqr_dl    += dx[i]*dx[i];
            }
            break;
        default:
            /* default to explicit simulation */
            return 1;
    }

    dl = sqrt(sqr_dl);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* molecule is coarse grained */
    if (dl > l2)
    {
        return 0;
    }
<<<<<<< HEAD
    else if (dl < adressr)
    {
        return 0;
=======
    /* molecule is explicit */
    else if (dl < adressr)
    {
        return 1;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    /* hybrid region */
    else
    {
<<<<<<< HEAD
        l=dl-adressr;
        H=adressw;
        H5=H*H*H*H*H;
        tmp=-30.0/H5*(l*l*l*l-2.0*l*l*l*H+l*l*H*H);
        if (dx[0]<0.0 && adresstype==eAdressXSplit)tmp*=-1.0;
        //printf ("DP %g %g\n", l, tmp);
        return tmp;
        //tmp=-2.0*cos((dl-adressr)*M_PI/2/adressw)*sin((dl-adressr)*M_PI/2/adressw)*M_PI/2/adressw;
        //return tmp;
=======
        tmp = cos((dl-adressr)*M_PI/2/adressw);
        return tmp*tmp;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}

void
update_adress_weights_com(FILE *               fplog,
                          int                  cg0,
                          int                  cg1,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc)
{
<<<<<<< HEAD
    int            icg,k,k0,k1,d;
    real           nrcg,inv_ncg,mtot,inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr,adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;
    real *         wfprime;
=======
    int            icg, k, k0, k1, d;
    real           nrcg, inv_ncg, mtot, inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr, adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    int n_hyb, n_ex, n_cg;

<<<<<<< HEAD
    n_hyb=0;
    n_cg=0;
    n_ex=0;
=======
    n_hyb = 0;
    n_cg  = 0;
    n_ex  = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
<<<<<<< HEAD
    wfprime                 = mdatoms->wfprime;
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    ref                = &(fr->adress_refs);


    /* Since this is center of mass AdResS, the vsite is not guaranteed
<<<<<<< HEAD
     * to be on the same node as the constructing atoms.  Therefore we 
=======
     * to be on the same node as the constructing atoms.  Therefore we
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
     * loop over the charge groups, calculate their center of mass,
     * then use this to calculate wf for each atom.  This wastes vsite
     * construction, but it's the only way to assure that the explicit
     * atoms have the same wf as their vsite. */

#ifdef DEBUG
<<<<<<< HEAD
    fprintf(fplog,"Calculating center of mass for charge groups %d to %d\n",
            cg0,cg1);
#endif
    cgindex = cgs->index;
    
    /* Compute the center of mass for all charge groups */
    for(icg=cg0; (icg<cg1); icg++) 
=======
    fprintf(fplog, "Calculating center of mass for charge groups %d to %d\n",
            cg0, cg1);
#endif
    cgindex = cgs->index;

    /* Compute the center of mass for all charge groups */
    for (icg = cg0; (icg < cg1); icg++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1-k0;
        if (nrcg == 1)
        {
<<<<<<< HEAD
            wf[k0] = adress_weight(x[k0],adresstype,adressr,adressw,ref,pbc,fr);
            wfprime[k0] = Dadress_weight(x[k0],adresstype,adressr,adressw,ref,pbc,fr);
            if (wf[k0]==0){ n_cg++;}
            else if (wf[k0]==1){ n_ex++;}
            else {n_hyb++;}

=======
            wf[k0] = adress_weight(x[k0], adresstype, adressr, adressw, ref, pbc, fr);
            if (wf[k0] == 0)
            {
                n_cg++;
            }
            else if (wf[k0] == 1)
            {
                n_ex++;
            }
            else
            {
                n_hyb++;
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        else
        {
            mtot = 0.0;
<<<<<<< HEAD
            for(k=k0; (k<k1); k++)
=======
            for (k = k0; (k < k1); k++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                mtot += massT[k];
            }
            if (mtot > 0.0)
            {
                inv_mtot = 1.0/mtot;
<<<<<<< HEAD
                
                clear_rvec(ix);
                for(k=k0; (k<k1); k++)
                {
                    for(d=0; (d<DIM); d++)
=======

                clear_rvec(ix);
                for (k = k0; (k < k1); k++)
                {
                    for (d = 0; (d < DIM); d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        ix[d] += x[k][d]*massT[k];
                    }
                }
<<<<<<< HEAD
                for(d=0; (d<DIM); d++)
=======
                for (d = 0; (d < DIM); d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    ix[d] *= inv_mtot;
                }
            }
            /* Calculate the center of gravity if the charge group mtot=0 (only vsites) */
            else
            {
                inv_ncg = 1.0/nrcg;

                clear_rvec(ix);
<<<<<<< HEAD
                for(k=k0; (k<k1); k++)
                {
                    for(d=0; (d<DIM); d++)
=======
                for (k = k0; (k < k1); k++)
                {
                    for (d = 0; (d < DIM); d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    {
                        ix[d] += x[k][d];
                    }
                }
<<<<<<< HEAD
                for(d=0; (d<DIM); d++)
=======
                for (d = 0; (d < DIM); d++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    ix[d] *= inv_ncg;
                }
            }

            /* Set wf of all atoms in charge group equal to wf of com */
<<<<<<< HEAD
            wf[k0] = adress_weight(ix,adresstype,adressr,adressw,ref,pbc, fr);
            wfprime[k0] = Dadress_weight(ix,adresstype,adressr,adressw,ref,pbc,fr);

            if (wf[k0]==0){ n_cg++;}
            else if (wf[k0]==1){ n_ex++;}
            else {n_hyb++;}

            for(k=(k0+1); (k<k1); k++)
            {
                wf[k] = wf[k0];
                wfprime[k]=wfprime[k0];
            }
        }
    }


    adress_set_kernel_flags(n_ex, n_hyb, n_cg, mdatoms);

    
}
void update_adress_weights_atom_per_atom(
                            int                  cg0,
                          int                  cg1,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc)
{
    int            icg,k,k0,k1,d;
    real           nrcg,inv_ncg,mtot,inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr,adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;
    real *         wfprime;
=======
            wf[k0] = adress_weight(ix, adresstype, adressr, adressw, ref, pbc, fr);

            if (wf[k0] == 0)
            {
                n_cg++;
            }
            else if (wf[k0] == 1)
            {
                n_ex++;
            }
            else
            {
                n_hyb++;
            }

            for (k = (k0+1); (k < k1); k++)
            {
                wf[k] = wf[k0];
            }
        }
    }
}

void update_adress_weights_atom_per_atom(
        int                  cg0,
        int                  cg1,
        t_block *            cgs,
        rvec                 x[],
        t_forcerec *         fr,
        t_mdatoms *          mdatoms,
        t_pbc *              pbc)
{
    int            icg, k, k0, k1, d;
    real           nrcg, inv_ncg, mtot, inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr, adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    int n_hyb, n_ex, n_cg;

<<<<<<< HEAD
    n_hyb=0;
    n_cg=0;
    n_ex=0;
=======
    n_hyb = 0;
    n_cg  = 0;
    n_ex  = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);
<<<<<<< HEAD
        wfprime              = mdatoms->wfprime;
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    cgindex = cgs->index;

    /* Weighting function is determined for each atom individually.
     * This is an approximation
     * as in the theory requires an interpolation based on the center of masses.
     * Should be used with caution */
<<<<<<< HEAD
   
    for (icg = cg0; (icg < cg1); icg++) {
        k0 = cgindex[icg];
        k1 = cgindex[icg + 1];
        nrcg = k1 - k0;

        for (k = (k0); (k < k1); k++) {
            wf[k] = adress_weight(x[k], adresstype, adressr, adressw, ref, pbc, fr);
            wfprime[k] = Dadress_weight(x[k],adresstype,adressr,adressw,ref,pbc,fr);
            if (wf[k] == 0) {
                n_cg++;
            } else if (wf[k] == 1) {
                n_ex++;
            } else {
=======

    for (icg = cg0; (icg < cg1); icg++)
    {
        k0   = cgindex[icg];
        k1   = cgindex[icg + 1];
        nrcg = k1 - k0;

        for (k = (k0); (k < k1); k++)
        {
            wf[k] = adress_weight(x[k], adresstype, adressr, adressw, ref, pbc, fr);
            if (wf[k] == 0)
            {
                n_cg++;
            }
            else if (wf[k] == 1)
            {
                n_ex++;
            }
            else
            {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                n_hyb++;
            }
        }

    }
<<<<<<< HEAD
    adress_set_kernel_flags(n_ex, n_hyb, n_cg, mdatoms);
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

void
update_adress_weights_cog(t_iparams            ip[],
                          t_ilist              ilist[],
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc)
{
<<<<<<< HEAD
    int            i,j,k,nr,nra,inc;
    int            ftype,adresstype;
    t_iatom        avsite,ai,aj,ak,al;
    t_iatom *      ia;
    real           adressr,adressw;
    rvec *         ref;
    real *         wf;
    int            n_hyb, n_ex, n_cg;
    
=======
    int            i, j, k, nr, nra, inc;
    int            ftype, adresstype;
    t_iatom        avsite, ai, aj, ak, al;
    t_iatom *      ia;
    real           adressr, adressw;
    rvec *         ref;
    real *         wf;
    int            n_hyb, n_ex, n_cg;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);


<<<<<<< HEAD
    n_hyb=0;
    n_cg=0;
    n_ex=0;
    
    gmx_fatal(FARGS,"adress_weights_cog not implemented for h-adress in %s, line %d",
                              __FILE__,__LINE__);
=======
    n_hyb = 0;
    n_cg  = 0;
    n_ex  = 0;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


    /* Since this is center of geometry AdResS, we know the vsite
     * is in the same charge group node as the constructing atoms.
     * Loop over vsite types, calculate the weight of the vsite,
     * then assign that weight to the constructing atoms. */

<<<<<<< HEAD
    for(ftype=0; (ftype<F_NRE); ftype++) 
    {
        if (interaction_function[ftype].flags & IF_VSITE) 
=======
    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            nra    = interaction_function[ftype].nratoms;
            nr     = ilist[ftype].nr;
            ia     = ilist[ftype].iatoms;
<<<<<<< HEAD
            
            for(i=0; (i<nr); ) 
=======

            for (i = 0; (i < nr); )
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            {
                /* The vsite and first constructing atom */
                avsite     = ia[1];
                ai         = ia[2];
<<<<<<< HEAD
                wf[avsite] = adress_weight(x[avsite],adresstype,adressr,adressw,ref,pbc,fr);
                wf[ai]     = wf[avsite];

                if (wf[ai]  == 0) {
                    n_cg++;
                } else if (wf[ai]  == 1) {
                    n_ex++;
                } else {
=======
                wf[avsite] = adress_weight(x[avsite], adresstype, adressr, adressw, ref, pbc, fr);
                wf[ai]     = wf[avsite];

                if (wf[ai]  == 0)
                {
                    n_cg++;
                }
                else if (wf[ai]  == 1)
                {
                    n_ex++;
                }
                else
                {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    n_hyb++;
                }

                /* Assign the vsite wf to rest of constructing atoms depending on type */
                inc = nra+1;
<<<<<<< HEAD
                switch (ftype) {
                case F_VSITE2:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    break;
                case F_VSITE3:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE3FD:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE3FAD:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE3OUT:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE4FD:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    al     = ia[5];
                    wf[al] = wf[avsite];
                    break;
                case F_VSITE4FDN:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    al     = ia[5];
                    wf[al] = wf[avsite];
                    break;
                case F_VSITEN:
                    inc    = 3*ip[ia[0]].vsiten.n;
                    for(j=3; j<inc; j+=3) 
                    {
                        ai = ia[j+2];
                        wf[ai] = wf[avsite];
                    }
                    break;
                default:
                    gmx_fatal(FARGS,"No such vsite type %d in %s, line %d",
                              ftype,__FILE__,__LINE__);
=======
                switch (ftype)
                {
                    case F_VSITE2:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        break;
                    case F_VSITE3:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE3FD:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE3FAD:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE3OUT:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE4FD:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        al     = ia[5];
                        wf[al] = wf[avsite];
                        break;
                    case F_VSITE4FDN:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        al     = ia[5];
                        wf[al] = wf[avsite];
                        break;
                    case F_VSITEN:
                        inc    = 3*ip[ia[0]].vsiten.n;
                        for (j = 3; j < inc; j += 3)
                        {
                            ai     = ia[j+2];
                            wf[ai] = wf[avsite];
                        }
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d",
                                  ftype, __FILE__, __LINE__);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
<<<<<<< HEAD

    adress_set_kernel_flags(n_ex, n_hyb, n_cg, mdatoms);
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

void
update_adress_weights_atom(int                  cg0,
                           int                  cg1,
                           t_block *            cgs,
                           rvec                 x[],
                           t_forcerec *         fr,
                           t_mdatoms *          mdatoms,
                           t_pbc *              pbc)
{
<<<<<<< HEAD
    int            icg,k,k0,k1;
    atom_id *      cgindex;
    int            adresstype;
    real           adressr,adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;   
    real *         wfprime;
    
=======
    int            icg, k, k0, k1;
    atom_id *      cgindex;
    int            adresstype;
    real           adressr, adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);
    cgindex            = cgs->index;
<<<<<<< HEAD
    wfprime              = mdatoms->wfprime;
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    /* Only use first atom in charge group.
     * We still can't be sure that the vsite and constructing
     * atoms are on the same processor, so we must calculate
     * in the same way as com adress. */
<<<<<<< HEAD
    
    for(icg=cg0; (icg<cg1); icg++) 
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        wf[k0] = adress_weight(x[k0],adresstype,adressr,adressw,ref,pbc,fr);
        wfprime[k0] = Dadress_weight(x[k0],adresstype,adressr,adressw,ref,pbc,fr);
        
        /* Set wf of all atoms in charge group equal to wf of first atom in charge group*/
        for(k=(k0+1); (k<k1); k++)
        {
            wf[k] = wf[k0];
            wfprime[k] =  wfprime[k0];
        }
    }
}

void adress_set_kernel_flags(int n_ex, int n_hyb, int n_cg, t_mdatoms * mdatoms){

    /* With domain decomposition we can check weather a cpu calculates only
     * coarse-grained or explicit interactions. If so we use standard gromacs kernels
     * on this proc. See also nonbonded.c */

    if (n_hyb ==0 && n_ex == 0){
     /* all particles on this proc are coarse-grained, use standard gromacs kernels */
        if (!mdatoms->purecg){
            mdatoms->purecg = TRUE;
           if (debug) fprintf (debug, "adress.c: pure cg kernels on this proc\n");
        }
    }
    else
    {
        if (mdatoms->purecg){
         /* now this processor has hybrid particles again, call the hybrid kernels */
            mdatoms->purecg = FALSE;
        }
    }

    if (n_hyb ==0 && n_cg == 0){
    /* all particles on this proc are atomistic, use standard gromacs kernels */
        if (!mdatoms->pureex){
             mdatoms->pureex = TRUE;
             if (debug) fprintf (debug, "adress.c: pure ex kernels on this proc\n");
        }
    }
    else
    {
        if (mdatoms->pureex){
            mdatoms->pureex = FALSE;
=======

    for (icg = cg0; (icg < cg1); icg++)
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        wf[k0]  = adress_weight(x[k0], adresstype, adressr, adressw, ref, pbc, fr);

        /* Set wf of all atoms in charge group equal to wf of first atom in charge group*/
        for (k = (k0+1); (k < k1); k++)
        {
            wf[k] = wf[k0];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
}

void
adress_thermo_force(int                  start,
                    int                  homenr,
                    t_block *            cgs,
                    rvec                 x[],
                    rvec                 f[],
                    t_forcerec *         fr,
                    t_mdatoms *          mdatoms,
<<<<<<< HEAD
                    t_pbc *              pbc, real * engdelta)
{
    int              iatom,n0,nnn,nrcg, i;
=======
                    t_pbc *              pbc)
{
    int              iatom, n0, nnn, nrcg, i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    int              adresstype;
    real             adressw, adressr;
    atom_id *        cgindex;
    unsigned short * ptype;
    rvec *           ref;
    real *           wf;
<<<<<<< HEAD
    real *           wfprime;
    real             tabscale;
    real *           ATFtab;
    rvec             dr;
    real             w,wsq,wmin1,wmin1sq,wp,wt,rinv, sqr_dl, dl;
    real             eps,eps2,F,Geps,Heps2,Fp,dmu_dwp,dwp_dr,fscal, VV, Y;
    

    adresstype       = fr->adress_type;
    adressw          = fr->adress_hy_width;
    adressr           = fr->adress_ex_width;
    cgindex          = cgs->index;
    ptype            = mdatoms->ptype;
    ref              = &(fr->adress_refs);
    wf               = mdatoms->wf;
    wfprime          = mdatoms->wfprime;
    *engdelta=0.0;
    

    for(iatom=start; (iatom<start+homenr); iatom++)
=======
    real             tabscale;
    real *           ATFtab;
    rvec             dr;
    real             w, wsq, wmin1, wmin1sq, wp, wt, rinv, sqr_dl, dl;
    real             eps, eps2, F, Geps, Heps2, Fp, dmu_dwp, dwp_dr, fscal;

    adresstype        = fr->adress_type;
    adressw           = fr->adress_hy_width;
    adressr           = fr->adress_ex_width;
    cgindex           = cgs->index;
    ptype             = mdatoms->ptype;
    ref               = &(fr->adress_refs);
    wf                = mdatoms->wf;

    for (iatom = start; (iatom < start+homenr); iatom++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        if (egp_coarsegrained(fr, mdatoms->cENER[iatom]))
        {
            if (ptype[iatom] == eptVSite)
            {
                w    = wf[iatom];
                /* is it hybrid or apply the thermodynamics force everywhere?*/
<<<<<<< HEAD
                if ( mdatoms->tf_table_index[iatom] != NO_TF_TABLE)
                {
                    if (fr->n_adress_tf_grps > 0 ){
                        /* multi component tf is on, select the right table */
                        ATFtab = fr->atf_tabs[mdatoms->tf_table_index[iatom]].tab;
                        tabscale = fr->atf_tabs[mdatoms->tf_table_index[iatom]].scale;
                        printf("atom %d table %d \n", iatom, mdatoms->tf_table_index[iatom]);
                    }
                    else {
                    /* just on component*/
                        ATFtab = fr->atf_tabs[DEFAULT_TF_TABLE].tab;
                        tabscale = fr->atf_tabs[DEFAULT_TF_TABLE].scale;
                        printf("atom %d NO COMPENSTATION \n", iatom);
                    }
                    
                    fscal            = 0;
                    sqr_dl =0.0;


                   if (pbc)
                    {
                        pbc_dx(pbc,(*ref),x[iatom],dr);
                    }
                    else
                    {
                        rvec_sub((*ref),x[iatom],dr);
                    }

                    switch(adresstype)
                    {
                    case eAdressXSplit:
                        /* plane through center of ref, varies in x direction */
                        sqr_dl         = dr[0]*dr[0];
                        rinv             = gmx_invsqrt(dr[0]*dr[0]);
                        break;
                    case eAdressSphere:
                        /* point at center of ref, assuming cubic geometry */
                        for(i=0;i<3;i++){
                            sqr_dl    += dr[i]*dr[i];
                        }
                        rinv             = gmx_invsqrt(sqr_dl);
                        break;
                    default:
                        /* This case should not happen */
                        rinv = 0.0;
                    }

                    
                    dl=wf[iatom];
=======
                if (mdatoms->tf_table_index[iatom] != NO_TF_TABLE)
                {
                    if (fr->n_adress_tf_grps > 0)
                    {
                        /* multi component tf is on, select the right table */
                        ATFtab   = fr->atf_tabs[mdatoms->tf_table_index[iatom]].data;
                        tabscale = fr->atf_tabs[mdatoms->tf_table_index[iatom]].scale;
                    }
                    else
                    {
                        /* just on component*/
                        ATFtab   = fr->atf_tabs[DEFAULT_TF_TABLE].data;
                        tabscale = fr->atf_tabs[DEFAULT_TF_TABLE].scale;
                    }

                    fscal            = 0;
                    if (pbc)
                    {
                        pbc_dx(pbc, (*ref), x[iatom], dr);
                    }
                    else
                    {
                        rvec_sub((*ref), x[iatom], dr);
                    }




                    /* calculate distace to adress center again */
                    sqr_dl = 0.0;
                    switch (adresstype)
                    {
                        case eAdressXSplit:
                            /* plane through center of ref, varies in x direction */
                            sqr_dl           = dr[0]*dr[0];
                            rinv             = gmx_invsqrt(dr[0]*dr[0]);
                            break;
                        case eAdressSphere:
                            /* point at center of ref, assuming cubic geometry */
                            for (i = 0; i < 3; i++)
                            {
                                sqr_dl    += dr[i]*dr[i];
                            }
                            rinv             = gmx_invsqrt(sqr_dl);
                            break;
                        default:
                            /* This case should not happen */
                            rinv = 0.0;
                    }

                    dl = sqrt(sqr_dl);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    /* table origin is adress center */
                    wt               = dl*tabscale;
                    n0               = wt;
                    eps              = wt-n0;
                    eps2             = eps*eps;
                    nnn              = 4*n0;
<<<<<<< HEAD
                    Y                = ATFtab[nnn];
=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    F                = ATFtab[nnn+1];
                    Geps             = eps*ATFtab[nnn+2];
                    Heps2            = eps2*ATFtab[nnn+3];
                    Fp               = F+Geps+Heps2;
                    F                = (Fp+Geps+2.0*Heps2)*tabscale;
<<<<<<< HEAD
                    VV               = Y+eps*Fp;

                    fscal            = F*wfprime[iatom];
                    *engdelta       += VV;

                    if (adresstype == eAdressXSplit){
                        f[iatom][0]        += fscal;
                    }
                    else if (adresstype == eAdressSphere)
                    {
                        f[iatom][0]    += fscal*dr[0]*rinv;
                        f[iatom][1]    += fscal*dr[1]*rinv;
                        f[iatom][2]    += fscal*dr[2]*rinv;
                    }

                }       
=======

                    fscal            = F*rinv;

                    f[iatom][0]        += fscal*dr[0];
                    if (adresstype != eAdressXSplit)
                    {
                        f[iatom][1]    += fscal*dr[1];
                        f[iatom][2]    += fscal*dr[2];
                    }
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
}

gmx_bool egp_explicit(t_forcerec *   fr, int egp_nr)
{
    return fr->adress_group_explicit[egp_nr];
}

gmx_bool egp_coarsegrained(t_forcerec *   fr, int egp_nr)
{
<<<<<<< HEAD
   return !fr->adress_group_explicit[egp_nr];
=======
    return !fr->adress_group_explicit[egp_nr];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}
