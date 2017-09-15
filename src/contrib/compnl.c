/*
<<<<<<< HEAD
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
 * GROwing Monsters And Cloning Shrimps
=======
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

#include "ns.h"
#include "smalloc.h"
#include "wnblist.h"
#include "futil.h"
#include "macros.h"
#include "statutil.h"
#include "copyrite.h"
#include "confio.h"
#include "pbc.h"
#include "vec.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]compnl[tt] compares two neighborlists as generated by [TT]mdrun[tt]",
    "in the log file, when the environment variable DUMPNL is set to",
    "a number larger than 0. [TT]compnl[tt] is mainly used for debugging the",
    "[TT]mdrun[tt] internals and not for end-users."
  };
  FILE    *in,*out;
  int     i,j,nmiss,mod;
  char    **fn,title[256];
  int     ***mat,nnb;
  real    mb;
  gmx_bool    bConf;
  rvec    *x = NULL;
  rvec    dx;
  matrix  box;
  t_atoms atoms;
  t_pbc   pbc;
  
  t_filenm fnm[] = {
    { efLOG, "-f1", NULL, ffREAD },
    { efLOG, "-f2", NULL, ffREAD },
    { efOUT, "-o",  "compnl", ffWRITE },
    { efSTX, "-c",  NULL, ffOPTRD }
  };
#define NFILE asize(fnm)
  static int natoms=648;
  static gmx_bool bSymm=TRUE;
  static t_pargs pa[] = {
    { "-nat",  FALSE, etINT, { &natoms }, "Number of atoms" },
    { "-symm", FALSE, etBOOL,{ &bSymm  }, "Symmetrize the matrices" },
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  bConf = (opt2bSet("-c",NFILE,fnm));
  if (bConf) {
    get_stx_coordnum (opt2fn("-c",NFILE,fnm),&natoms);
    init_t_atoms(&atoms,natoms,FALSE);
    snew(x,natoms);
    read_stx_conf(opt2fn("-c",NFILE,fnm),title,&atoms,x,NULL,box);
    set_pbc(&pbc,box);
  }
  snew(fn,2);
  fn[0] = opt2fn("-f1",NFILE,fnm);
  fn[1] = opt2fn("-f2",NFILE,fnm);
  
  snew(mat,2);  
  out = gmx_fio_fopen(ftp2fn(efOUT,NFILE,fnm),"w");
  mb  = sizeof(int)*sqr(natoms/1024.0);
  for(i=0; (i<2); i++) {
    in = gmx_fio_fopen(fn[i],"r");
    fprintf(stderr,"Reading %s\n",fn[i]);
    fprintf(out,   "Reading %s\n",fn[i]);
    fprintf(stderr,"Going to allocate %.0f Mb of memory\n",mb);
    fprintf(out,   "Going to allocate %.0f Mb of memory\n",mb);
    snew(mat[i],natoms);
    for(j=0; (j<natoms); j++) 
      snew(mat[i][j],natoms);
    nnb = read_nblist(in,out,mat[i],natoms,bSymm);
    gmx_fio_fclose(in);
    fprintf(stderr,"Interaction matrix %d has %d entries\n",i,nnb);
    fprintf(out,   "Interaction matrix %d has %d entries\n",i,nnb);
  }
  fprintf(stderr,"Comparing Interaction Matrices\n");
  mod=1;
  nmiss = 0;
  for(i=0; (i<natoms); i+=mod) {
    for(j=0; (j<natoms); j+=mod) {
      if (mat[0][i][j] != mat[1][i][j]) {
	fprintf(out,"i: %5d, j: %5d, shift[%s]: %3d, shift[%s]: %3d",
		i,j,fn[0],mat[0][i][j]-1,fn[1],mat[1][i][j]-1);
	if (bConf) {
	  pbc_dx(&pbc,x[i],x[j],dx);
	  fprintf(out," dist: %8.3f\n",norm(dx));
	}
	else
	  fprintf(out,"\n");
	nmiss++;
      }
    }
  }
  fprintf(out,"There were %d mismatches\n",nmiss);
  fprintf(out,"Done.\n");
  gmx_fio_fclose(out);
  fprintf(stderr,"There were %d mismatches\n",nmiss);
  fprintf(stderr,"Finished\n");
  
  thanx(stdout);
  
  return 0;
}



