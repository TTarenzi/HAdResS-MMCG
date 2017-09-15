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
 * GROningen Mixture of Alchemy and Childrens' Stories
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

#include <string.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "filenm.h"
#include "futil.h"
#include "wman.h"
#include "xdrf.h"
#include "macros.h"

<<<<<<< HEAD
#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

/* NOTE: this was a cesspool of thread-unsafe code, has now been 
 properly proteced by mutexes (hopefully). */

/* XDR should be available on all platforms now, 
=======
#ifdef GMX_THREAD_MPI
#include "thread_mpi.h"
#endif

/* NOTE: this was a cesspool of thread-unsafe code, has now been
   properly proteced by mutexes (hopefully). */

/* XDR should be available on all platforms now,
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * but we keep the possibility of turning it off...
 */
#define USE_XDR

/* Use bitflag ... */
#define IS_SET(fn) ((fn.flag & ffSET) != 0)
#define IS_OPT(fn) ((fn.flag & ffOPT) != 0)
#define IS_MULT(fn) ((fn.flag & ffMULT) != 0)
#define UN_SET(fn) (fn.flag = (fn.flag & ~ffSET))
#define DO_SET(fn) (fn.flag = (fn.flag |  ffSET))

enum
{
    eftASC, eftBIN, eftXDR, eftGEN, eftNR
};

<<<<<<< HEAD
/* To support multiple file types with one general (eg TRX) we have 
 * these arrays.
 */
static const int trxs[] =
    {
#ifdef USE_XDR 
        efXTC, efTRR, efCPT,
#endif
        efTRJ, efGRO, efG96, efPDB, efG87 };
#define NTRXS asize(trxs)

static const int tros[] =
    {
#ifdef USE_XDR 
        efXTC, efTRR,
#endif
        efTRJ, efGRO, efG96, efPDB, efG87 };
#define NTROS asize(tros)

static const int trns[] =
    {
#ifdef USE_XDR
        efTRR, efCPT,
#endif
        efTRJ };
#define NTRNS asize(trns)

static const int stos[] =
    { efGRO, efG96, efPDB, efBRK, efENT, efESP, efXYZ };
#define NSTOS asize(stos)

static const int stxs[] =
    { efGRO, efG96, efPDB, efBRK, efENT, efESP, efXYZ,
#ifdef USE_XDR 
        efTPR,
#endif 
        efTPB, efTPA };
#define NSTXS asize(stxs)

static const int tpxs[] =
    {
#ifdef USE_XDR
        efTPR,
#endif
        efTPB, efTPA };
#define NTPXS asize(tpxs)

static const int tpss[] =
    {
#ifdef USE_XDR
        efTPR,
#endif
        efTPB, efTPA, efGRO, efG96, efPDB, efBRK, efENT };
=======
/* To support multiple file types with one general (eg TRX) we have
 * these arrays.
 */
static const int trxs[] =
{
#ifdef USE_XDR
    efXTC, efTRR, efCPT,
#endif
    efTRJ, efGRO, efG96, efPDB, efG87
};
#define NTRXS asize(trxs)

static const int tros[] =
{
#ifdef USE_XDR
    efXTC, efTRR,
#endif
    efTRJ, efGRO, efG96, efPDB, efG87
};
#define NTROS asize(tros)

static const int trns[] =
{
#ifdef USE_XDR
    efTRR, efCPT,
#endif
    efTRJ
};
#define NTRNS asize(trns)

static const int stos[] =
{ efGRO, efG96, efPDB, efBRK, efENT, efESP, efXYZ };
#define NSTOS asize(stos)

static const int stxs[] =
{
    efGRO, efG96, efPDB, efBRK, efENT, efESP, efXYZ,
#ifdef USE_XDR
    efTPR,
#endif
    efTPB, efTPA
};
#define NSTXS asize(stxs)

static const int tpxs[] =
{
#ifdef USE_XDR
    efTPR,
#endif
    efTPB, efTPA
};
#define NTPXS asize(tpxs)

static const int tpss[] =
{
#ifdef USE_XDR
    efTPR,
#endif
    efTPB, efTPA, efGRO, efG96, efPDB, efBRK, efENT
};
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#define NTPSS asize(tpss)

typedef struct
{
<<<<<<< HEAD
    int ftype;
=======
    int         ftype;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    const char *ext;
    const char *defnm;
    const char *defopt;
    const char *descr;
<<<<<<< HEAD
    int ntps;
    const int *tps;
=======
    int         ntps;
    const int  *tps;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
} t_deffile;

/* this array should correspond to the enum in include/types/filenm.h */
static const t_deffile
    deffile[efNR] =
{
    { eftASC, ".mdp", "grompp", "-f", "grompp input file with MD parameters" },
    { eftASC, ".gct", "gct",    "-f", "General coupling stuff"},
    { eftGEN, ".???", "traj", "-f",
      "Trajectory: xtc trr trj gro g96 pdb cpt", NTRXS, trxs },
    { eftGEN, ".???", "trajout", "-f",
      "Trajectory: xtc trr trj gro g96 pdb", NTROS, tros },
    { eftGEN, ".???", "traj", NULL,
      "Full precision trajectory: trr trj cpt", NTRNS, trns },
    { eftXDR, ".trr", "traj", NULL, "Trajectory in portable xdr format" },
    { eftBIN, ".trj", "traj", NULL, "Trajectory file (architecture specific)" },
    { eftXDR, ".xtc", "traj", NULL,
      "Compressed trajectory (portable xdr format)" },
    { eftASC, ".g87", "gtraj", NULL, "Gromos-87 ASCII trajectory format" },
    { eftXDR, ".edr", "ener",   NULL, "Energy file"},
<<<<<<< HEAD
    { eftGEN, ".???", "conf", "-c", "Structure file: gro g96 pdb tpr etc.", 
      NSTXS, stxs },
    { eftGEN, ".???", "out", "-o", "Structure file: gro g96 pdb etc.", 
=======
    { eftGEN, ".???", "conf", "-c", "Structure file: gro g96 pdb tpr etc.",
      NSTXS, stxs },
    { eftGEN, ".???", "out", "-o", "Structure file: gro g96 pdb etc.",
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
      NSTOS, stos },
    { eftASC, ".gro", "conf", "-c", "Coordinate file in Gromos-87 format" },
    { eftASC, ".g96", "conf", "-c", "Coordinate file in Gromos-96 format" },
    { eftASC, ".pdb", "eiwit",  "-f", "Protein data bank file"},
    { eftASC, ".brk", "eiwit",  "-f", "Brookhaven data bank file"},
    { eftASC, ".ent", "eiwit", "-f", "Entry in the protein date bank" },
    { eftASC, ".esp", "conf", "-f", "Coordinate file in Espresso format" },
    { eftASC, ".pqr", "state",  "-o", "Coordinate file for MEAD"},
    { eftASC, ".xyz", "conf", "-o", "Coordinate file for some other programs" },
<<<<<<< HEAD
    { eftXDR, ".cpt", "state",  "-cp","Checkpoint file"},
    { eftASC, ".log", "run",    "-l", "Log file"},
    { eftASC, ".xvg", "graph",  "-o", "xvgr/xmgr file"},
    { eftASC, ".out", "hello",  "-o", "Generic output file"},
    { eftASC, ".ndx", "index",  "-n", "Index file",},
=======
    { eftXDR, ".cpt", "state",  "-cp", "Checkpoint file"},
    { eftASC, ".log", "run",    "-l", "Log file"},
    { eftASC, ".xvg", "graph",  "-o", "xvgr/xmgr file"},
    { eftASC, ".out", "hello",  "-o", "Generic output file"},
    { eftASC, ".ndx", "index",  "-n", "Index file", },
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    { eftASC, ".top", "topol",  "-p", "Topology file"},
    { eftASC, ".itp", "topinc", NULL, "Include file for topology"},
    { eftGEN, ".???", "topol", "-s", "Run input file: tpr tpb tpa",
      NTPXS, tpxs },
    { eftGEN, ".???", "topol", "-s",
      "Structure+mass(db): tpr tpb tpa gro g96 pdb", NTPSS, tpss },
    { eftXDR, ".tpr", "topol",  "-s", "Portable xdr run input file"},
    { eftASC, ".tpa", "topol",  "-s", "Ascii run input file"},
    { eftBIN, ".tpb", "topol",  "-s", "Binary run input file"},
    { eftASC, ".tex", "doc",    "-o", "LaTeX file"},
    { eftASC, ".rtp", "residue", NULL, "Residue Type file used by pdb2gmx" },
    { eftASC, ".atp", "atomtp", NULL, "Atomtype file used by pdb2gmx" },
    { eftASC, ".hdb", "polar",  NULL, "Hydrogen data base"},
    { eftASC, ".dat", "nnnice", NULL, "Generic data file"},
    { eftASC, ".dlg", "user",   NULL, "Dialog Box data for ngmx"},
    { eftASC, ".map", "ss", NULL, "File that maps matrix data to colors" },
    { eftASC, ".eps", "plot", NULL, "Encapsulated PostScript (tm) file" },
    { eftASC, ".mat", "ss",     NULL, "Matrix Data file"},
    { eftASC, ".m2p", "ps",     NULL, "Input file for mat2ps"},
<<<<<<< HEAD
    { eftXDR, ".mtx", "hessian","-m", "Hessian matrix"},
    { eftASC, ".edi", "sam",    NULL, "ED sampling input"},
    { eftASC, ".edo", "sam",    NULL, "ED sampling output"},
    { eftASC, ".hat", "gk", NULL, "Fourier transform of spread function" },
    { eftASC, ".cub", "pot",  NULL, "Gaussian cube file" },
    { eftASC, ".xpm", "root", NULL, "X PixMap compatible matrix file" },
    { eftASC, "", "rundir", NULL, "Run directory" } 
=======
    { eftXDR, ".mtx", "hessian", "-m", "Hessian matrix"},
    { eftASC, ".edi", "sam",    NULL, "ED sampling input"},
    { eftASC, ".hat", "gk", NULL, "Fourier transform of spread function" },
    { eftASC, ".cub", "pot",  NULL, "Gaussian cube file" },
    { eftASC, ".xpm", "root", NULL, "X PixMap compatible matrix file" },
    { eftASC, "", "rundir", NULL, "Run directory" }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
};

static char *default_file_name = NULL;

<<<<<<< HEAD
#ifdef GMX_THREADS
static tMPI_Thread_mutex_t filenm_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
=======
#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t filenm_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif

#define NZEXT 2
const char *z_ext[NZEXT] =
<<<<<<< HEAD
    { ".gz", ".Z" };
=======
{ ".gz", ".Z" };
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

void set_default_file_name(const char *name)
{
    int i;
<<<<<<< HEAD
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&filenm_mutex);
#endif
    default_file_name = strdup(name);
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&filenm_mutex);
#endif
    default_file_name = strdup(name);
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&filenm_mutex);
#endif

#if 0
<<<<<<< HEAD
    for(i=0; i<efNR; i++)
    deffile[i].defnm = default_file_name;
=======
    for (i = 0; i < efNR; i++)
    {
        deffile[i].defnm = default_file_name;
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#endif
}

const char *ftp2ext(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
<<<<<<< HEAD
        return deffile[ftp].ext + 1;
    else
        return "unknown";
=======
    {
        return deffile[ftp].ext + 1;
    }
    else
    {
        return "unknown";
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

const char *ftp2ext_generic(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        switch (ftp)
        {
<<<<<<< HEAD
        case efTRX:
            return "trx";
        case efTRN:
            return "trn";
        case efSTO:
            return "sto";
        case efSTX:
            return "stx";
        case efTPX:
            return "tpx";
        case efTPS:
            return "tps";
        default:
            return ftp2ext(ftp);
        }
    }
    else
        return "unknown";
=======
            case efTRX:
                return "trx";
            case efTRN:
                return "trn";
            case efSTO:
                return "sto";
            case efSTX:
                return "stx";
            case efTPX:
                return "tpx";
            case efTPS:
                return "tps";
            default:
                return ftp2ext(ftp);
        }
    }
    else
    {
        return "unknown";
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

const char *ftp2desc(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
<<<<<<< HEAD
        return deffile[ftp].descr;
    else
        return "unknown filetype";
=======
    {
        return deffile[ftp].descr;
    }
    else
    {
        return "unknown filetype";
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

const char *ftp2ftype(int ftp)
{
    if ((ftp >= 0) && (ftp < efNR))
    {
        switch (deffile[ftp].ftype)
        {
<<<<<<< HEAD
        case eftASC:
            return "ASCII";
        case eftBIN:
            return "Binary";
        case eftXDR:
            return "XDR portable";
        case eftGEN:
            return "";
        default:
            gmx_fatal(FARGS, "Unknown filetype %d in ftp2ftype",deffile[ftp].ftype);
            break;
=======
            case eftASC:
                return "ASCII";
            case eftBIN:
                return "Binary";
            case eftXDR:
                return "XDR portable";
            case eftGEN:
                return "";
            default:
                gmx_fatal(FARGS, "Unknown filetype %d in ftp2ftype", deffile[ftp].ftype);
                break;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
    }
    return "unknown";
}

const char *ftp2defnm(int ftp)
{
    const char *buf = NULL;

<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_lock(&filenm_mutex);
#endif

    if (default_file_name)
    {
        buf = default_file_name;
    }
    else
    {
        if ((0 <= ftp) && (ftp < efNR))
        {
            buf = deffile[ftp].defnm;
        }
    }
<<<<<<< HEAD
#ifdef GMX_THREADS
=======
#ifdef GMX_THREAD_MPI
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    tMPI_Thread_mutex_unlock(&filenm_mutex);
#endif

    return buf;
}

void pr_def(FILE *fp, int ftp)
{
    const t_deffile *df;
<<<<<<< HEAD
    const char *s = NULL;
    char *flst, *tmp, *desc;
    const char *ext;
    const char *defnm;

    df = &(deffile[ftp]);
=======
    const char      *s = NULL;
    char            *flst, *tmp, *desc;
    const char      *ext;
    const char      *defnm;

    df    = &(deffile[ftp]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    defnm = ftp2defnm(ftp);
    /* find default file extension and \tt-ify description */
    /* FIXME: The constness should not be cast away */
    flst = (char *) "";
    desc = strdup(df->descr);

    if (df->ntps)
    {
        ext = deffile[df->tps[0]].ext;
        tmp = strstr(desc, ": ") + 1;
        if (tmp)
        {
            tmp[0] = '\0';
            tmp++;
<<<<<<< HEAD
            snew(flst,strlen(tmp)+6);
=======
            snew(flst, strlen(tmp)+6);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            strcpy(flst, " \\tt ");
            strcat(flst, tmp);
        }
    }
    else
    {
        ext = df->ext;
    }
    /* now skip dot */
    if (ext[0])
<<<<<<< HEAD
        ext++;
    else
        ext = "";
    /* set file contents type */
    switch (df->ftype)
    {
    case eftASC:
        s = "Asc";
        break;
    case eftBIN:
        s = "Bin";
        break;
    case eftXDR:
        s = "xdr";
        break;
    case eftGEN:
        s = "";
        break;
    default:
        gmx_fatal(FARGS, "Unimplemented filetype %d %d",ftp,
        df->ftype);
    }
    fprintf(fp,"\\tt %8s & \\tt %3s & %3s & \\tt %2s & %s%s \\\\[-0.1ex]\n",
        defnm, ext, s, df->defopt ? df->defopt : "",
        check_tex(desc),check_tex(flst));
=======
    {
        ext++;
    }
    else
    {
        ext = "";
    }
    /* set file contents type */
    switch (df->ftype)
    {
        case eftASC:
            s = "Asc";
            break;
        case eftBIN:
            s = "Bin";
            break;
        case eftXDR:
            s = "xdr";
            break;
        case eftGEN:
            s = "";
            break;
        default:
            gmx_fatal(FARGS, "Unimplemented filetype %d %d", ftp,
                      df->ftype);
    }
    fprintf(fp, "\\tt %8s & \\tt %3s & %3s & \\tt %2s & %s%s \\\\[-0.1ex]\n",
            defnm, ext, s, df->defopt ? df->defopt : "",
            check_tex(desc), check_tex(flst));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    free(desc);
}

void pr_fns(FILE *fp, int nf, const t_filenm tfn[])
{
<<<<<<< HEAD
    int i, f;
    size_t j;
    char buf[256], *wbuf, opt_buf[32];
=======
    int    i, f;
    size_t j;
    char   buf[256], *wbuf, opt_buf[32];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#define OPTLEN 4
#define NAMELEN 14
    fprintf(fp, "%6s %12s  %-12s %s\n", "Option", "Filename", "Type",
            "Description");
    fprintf(fp,
            "------------------------------------------------------------\n");
    for (i = 0; (i < nf); i++)
    {
        for (f = 0; (f < tfn[i].nfiles); f++)
        {
            sprintf(buf, "%4s %14s  %-12s ", (f == 0) ? tfn[i].opt : "",
                    tfn[i].fns[f], (f == 0) ? fileopt(tfn[i].flag, opt_buf, 32)
<<<<<<< HEAD
                        : "");
            if (f < tfn[i].nfiles - 1)
                fprintf(fp, "%s\n", buf);
=======
                    : "");
            if (f < tfn[i].nfiles - 1)
            {
                fprintf(fp, "%s\n", buf);
            }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        }
        if (tfn[i].nfiles > 0)
        {
            strcat(buf, deffile[tfn[i].ftp].descr);
            if ((strlen(tfn[i].opt) > OPTLEN)
                && (strlen(tfn[i].opt) <= ((OPTLEN + NAMELEN)
<<<<<<< HEAD
                    - strlen(tfn[i].fns[tfn[i].nfiles - 1]))))
            {
                for (j = strlen(tfn[i].opt); j < strlen(buf)
                    - (strlen(tfn[i].opt) - OPTLEN) + 1; j++)
                    buf[j] = buf[j + strlen(tfn[i].opt) - OPTLEN];
=======
                                           - strlen(tfn[i].fns[tfn[i].nfiles - 1]))))
            {
                for (j = strlen(tfn[i].opt); j < strlen(buf)
                     - (strlen(tfn[i].opt) - OPTLEN) + 1; j++)
                {
                    buf[j] = buf[j + strlen(tfn[i].opt) - OPTLEN];
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            wbuf = wrap_lines(buf, 78, 35, FALSE);
            fprintf(fp, "%s\n", wbuf);
            sfree(wbuf);
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
}

void pr_fopts(FILE *fp, int nf, const t_filenm tfn[], int shell)
{
    int i, j;

    switch (shell)
    {
<<<<<<< HEAD
    case eshellCSH:
        for (i = 0; (i < nf); i++)
        {
            fprintf(fp, " \"n/%s/f:*.", tfn[i].opt);
            if (deffile[tfn[i].ftp].ntps)
            {
                fprintf(fp, "{");
                for (j = 0; j < deffile[tfn[i].ftp].ntps; j++)
                {
                    if (j > 0)
                        fprintf(fp, ",");
                    fprintf(fp, "%s", deffile[deffile[tfn[i].ftp].tps[j]].ext
                            + 1);
                }
                fprintf(fp, "}");
            }
            else
                fprintf(fp, "%s", deffile[tfn[i].ftp].ext + 1);
            fprintf(fp, "{");
            for (j = 0; j < NZEXT; j++)
                fprintf(fp, ",%s", z_ext[j]);
            fprintf(fp, "}/\"");
        }
        break;
    case eshellBASH:
        for (i = 0; (i < nf); i++)
        {
            fprintf(fp, "%s) COMPREPLY=( $(compgen -X '!*.", tfn[i].opt);
            if (deffile[tfn[i].ftp].ntps)
            {
                fprintf(fp, "+(");
                for (j = 0; j < deffile[tfn[i].ftp].ntps; j++)
                {
                    if (j > 0)
                        fprintf(fp, "|");
                    fprintf(fp, "%s", deffile[deffile[tfn[i].ftp].tps[j]].ext
                            + 1);
                }
                fprintf(fp, ")");
            }
            else
                fprintf(fp, "%s", deffile[tfn[i].ftp].ext + 1);
            fprintf(fp, "*(");
            for (j = 0; j < NZEXT; j++)
            {
                if (j > 0)
                    fprintf(fp, "|");
                fprintf(fp, "%s", z_ext[j]);
            }
            fprintf(fp, ")' -f $c ; compgen -S '/' -X '.*' -d $c ));;\n");
        }
        break;
    case eshellZSH:
        for (i = 0; (i < nf); i++)
        {
            fprintf(fp, "- 'c[-1,%s]' -g '*.", tfn[i].opt);
            if (deffile[tfn[i].ftp].ntps)
            {
                fprintf(fp, "(");
                for (j = 0; j < deffile[tfn[i].ftp].ntps; j++)
                {
                    if (j > 0)
                        fprintf(fp, "|");
                    fprintf(fp, "%s", deffile[deffile[tfn[i].ftp].tps[j]].ext
                            + 1);
                }
                fprintf(fp, ")");
            }
            else
                fprintf(fp, "%s", deffile[tfn[i].ftp].ext + 1);
            fprintf(fp, "(");
            for (j = 0; j < NZEXT; j++)
                fprintf(fp, "|%s", z_ext[j]);
            fprintf(fp, ") *(/)' ");
        }
        break;
=======
        case eshellCSH:
            for (i = 0; (i < nf); i++)
            {
                fprintf(fp, " \"n/%s/f:*.", tfn[i].opt);
                if (deffile[tfn[i].ftp].ntps)
                {
                    fprintf(fp, "{");
                    for (j = 0; j < deffile[tfn[i].ftp].ntps; j++)
                    {
                        if (j > 0)
                        {
                            fprintf(fp, ",");
                        }
                        fprintf(fp, "%s", deffile[deffile[tfn[i].ftp].tps[j]].ext
                                + 1);
                    }
                    fprintf(fp, "}");
                }
                else
                {
                    fprintf(fp, "%s", deffile[tfn[i].ftp].ext + 1);
                }
                fprintf(fp, "{");
                for (j = 0; j < NZEXT; j++)
                {
                    fprintf(fp, ",%s", z_ext[j]);
                }
                fprintf(fp, "}/\"");
            }
            break;
        case eshellBASH:
            for (i = 0; (i < nf); i++)
            {
                fprintf(fp, "%s) COMPREPLY=( $(compgen -X '!*.", tfn[i].opt);
                if (deffile[tfn[i].ftp].ntps)
                {
                    fprintf(fp, "+(");
                    for (j = 0; j < deffile[tfn[i].ftp].ntps; j++)
                    {
                        if (j > 0)
                        {
                            fprintf(fp, "|");
                        }
                        fprintf(fp, "%s", deffile[deffile[tfn[i].ftp].tps[j]].ext
                                + 1);
                    }
                    fprintf(fp, ")");
                }
                else
                {
                    fprintf(fp, "%s", deffile[tfn[i].ftp].ext + 1);
                }
                fprintf(fp, "*(");
                for (j = 0; j < NZEXT; j++)
                {
                    if (j > 0)
                    {
                        fprintf(fp, "|");
                    }
                    fprintf(fp, "%s", z_ext[j]);
                }
                fprintf(fp, ")' -f $c ; compgen -S '/' -X '.*' -d $c ));;\n");
            }
            break;
        case eshellZSH:
            for (i = 0; (i < nf); i++)
            {
                fprintf(fp, "- 'c[-1,%s]' -g '*.", tfn[i].opt);
                if (deffile[tfn[i].ftp].ntps)
                {
                    fprintf(fp, "(");
                    for (j = 0; j < deffile[tfn[i].ftp].ntps; j++)
                    {
                        if (j > 0)
                        {
                            fprintf(fp, "|");
                        }
                        fprintf(fp, "%s", deffile[deffile[tfn[i].ftp].tps[j]].ext
                                + 1);
                    }
                    fprintf(fp, ")");
                }
                else
                {
                    fprintf(fp, "%s", deffile[tfn[i].ftp].ext + 1);
                }
                fprintf(fp, "(");
                for (j = 0; j < NZEXT; j++)
                {
                    fprintf(fp, "|%s", z_ext[j]);
                }
                fprintf(fp, ") *(/)' ");
            }
            break;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
}

static void check_opts(int nf, t_filenm fnm[])
{
<<<<<<< HEAD
    int i;
=======
    int              i;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    const t_deffile *df;

    for (i = 0; (i < nf); i++)
    {
        df = &(deffile[fnm[i].ftp]);
        if (fnm[i].opt == NULL)
        {
            if (df->defopt == NULL)
            {
                gmx_fatal(FARGS, "No default cmd-line option for %s (type %d)\n",
<<<<<<< HEAD
                          deffile[fnm[i].ftp].ext,fnm[i].ftp);
            }
            else
            {
                fnm[i].opt=df->defopt;
=======
                          deffile[fnm[i].ftp].ext, fnm[i].ftp);
            }
            else
            {
                fnm[i].opt = df->defopt;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }
}

int fn2ftp(const char *fn)
{
<<<<<<< HEAD
    int i, len;
=======
    int         i, len;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    const char *feptr;
    const char *eptr;

    if (!fn)
<<<<<<< HEAD
        return efNR;

    len = strlen(fn);
    if ((len >= 4) && (fn[len - 4] == '.'))
        feptr = &(fn[len - 4]);
    else
        return efNR;

    for (i = 0; (i < efNR); i++)
        if ((eptr = deffile[i].ext) != NULL)
            if (gmx_strcasecmp(feptr, eptr) == 0)
                break;
=======
    {
        return efNR;
    }

    len = strlen(fn);
    if ((len >= 4) && (fn[len - 4] == '.'))
    {
        feptr = &(fn[len - 4]);
    }
    else
    {
        return efNR;
    }

    for (i = 0; (i < efNR); i++)
    {
        if ((eptr = deffile[i].ext) != NULL)
        {
            if (gmx_strcasecmp(feptr, eptr) == 0)
            {
                break;
            }
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    return i;
}

static void set_extension(char *buf, int ftp)
{
<<<<<<< HEAD
    int len, extlen;
    const t_deffile *df;

    /* check if extension is already at end of filename */
    df = &(deffile[ftp]);
    len = strlen(buf);
    extlen = strlen(df->ext);
    if ((len <= extlen) || (gmx_strcasecmp(&(buf[len - extlen]), df->ext) != 0))
        strcat(buf, df->ext);
=======
    int              len, extlen;
    const t_deffile *df;

    /* check if extension is already at end of filename */
    df     = &(deffile[ftp]);
    len    = strlen(buf);
    extlen = strlen(df->ext);
    if ((len <= extlen) || (gmx_strcasecmp(&(buf[len - extlen]), df->ext) != 0))
    {
        strcat(buf, df->ext);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

static void add_filenm(t_filenm *fnm, const char *filenm)
{
    srenew(fnm->fns, fnm->nfiles+1);
    fnm->fns[fnm->nfiles] = strdup(filenm);
    fnm->nfiles++;
}

static void set_grpfnm(t_filenm *fnm, const char *name, gmx_bool bCanNotOverride)
{
<<<<<<< HEAD
    char buf[256], buf2[256];
    int i, type;
    gmx_bool bValidExt;
    int nopts;
    const int *ftps;

    nopts = deffile[fnm->ftp].ntps;
    ftps = deffile[fnm->ftp].tps;
    if ((nopts == 0) || (ftps == NULL))
        gmx_fatal(FARGS, "nopts == 0 || ftps == NULL");
=======
    char       buf[256], buf2[256];
    int        i, type;
    gmx_bool   bValidExt;
    int        nopts;
    const int *ftps;

    nopts = deffile[fnm->ftp].ntps;
    ftps  = deffile[fnm->ftp].tps;
    if ((nopts == 0) || (ftps == NULL))
    {
        gmx_fatal(FARGS, "nopts == 0 || ftps == NULL");
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    bValidExt = FALSE;
    if (name && (bCanNotOverride || (default_file_name == NULL)))
    {
<<<<<<< HEAD
        strcpy(buf,name);
=======
        strcpy(buf, name);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        /* First check whether we have a valid filename already */
        type = fn2ftp(name);
        if ((fnm->flag & ffREAD) && (fnm->ftp == efTRX))
        {
            /*if file exist don't add an extension for trajectory reading*/
<<<<<<< HEAD
            bValidExt = gmx_fexist(name); 
        }
        for(i=0; (i<nopts) && !bValidExt; i++)
=======
            bValidExt = gmx_fexist(name);
        }
        for (i = 0; (i < nopts) && !bValidExt; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            if (type == ftps[i])
            {
                bValidExt = TRUE;
            }
        }
    }
    else
<<<<<<< HEAD
        /* No name given, set the default name */
        strcpy(buf,ftp2defnm(fnm->ftp));
=======
    {
        /* No name given, set the default name */
        strcpy(buf, ftp2defnm(fnm->ftp));
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if (!bValidExt && (fnm->flag & ffREAD))
    {
        /* for input-files only: search for filenames in the directory */
<<<<<<< HEAD
        for(i=0; (i<nopts) && !bValidExt; i++)
        {
            type = ftps[i];
            strcpy(buf2,buf);
            set_extension(buf2,type);
            if (gmx_fexist(buf2))
            {
                bValidExt = TRUE;
                strcpy(buf,buf2);
=======
        for (i = 0; (i < nopts) && !bValidExt; i++)
        {
            type = ftps[i];
            strcpy(buf2, buf);
            set_extension(buf2, type);
            if (gmx_fexist(buf2))
            {
                bValidExt = TRUE;
                strcpy(buf, buf2);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
        }
    }

    if (!bValidExt)
    {
        /* Use the first extension type */
<<<<<<< HEAD
        set_extension(buf,ftps[0]);
=======
        set_extension(buf, ftps[0]);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }

    add_filenm(fnm, buf);
}

static void set_filenm(t_filenm *fnm, const char *name, gmx_bool bCanNotOverride,
                       gmx_bool bReadNode)
{
<<<<<<< HEAD
    /* Set the default filename, extension and option for those fields that 
=======
    /* Set the default filename, extension and option for those fields that
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
     * are not already set. An extension is added if not present, if fn = NULL
     * or empty, the default filename is given.
     */
    char buf[256];
<<<<<<< HEAD
    int i, len, extlen;
=======
    int  i, len, extlen;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    if ((fnm->flag & ffREAD) && !bReadNode)
    {
        return;
    }

    if ((fnm->ftp < 0) || (fnm->ftp >= efNR))
<<<<<<< HEAD
        gmx_fatal(FARGS, "file type out of range (%d)",fnm->ftp);

    if (name)
        strcpy(buf, name);
=======
    {
        gmx_fatal(FARGS, "file type out of range (%d)", fnm->ftp);
    }

    if (name)
    {
        strcpy(buf, name);
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    if ((fnm->flag & ffREAD) && name && gmx_fexist(name))
    {
        /* check if filename ends in .gz or .Z, if so remove that: */
        len = strlen(name);
<<<<<<< HEAD
        for (i=0; i<NZEXT; i++)
=======
        for (i = 0; i < NZEXT; i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            extlen = strlen(z_ext[i]);
            if (len > extlen)
            {
<<<<<<< HEAD
                if (gmx_strcasecmp(name+len-extlen,z_ext[i]) == 0)
                {
                    buf[len-extlen]='\0';
=======
                if (gmx_strcasecmp(name+len-extlen, z_ext[i]) == 0)
                {
                    buf[len-extlen] = '\0';
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    break;
                }
            }
        }
    }

    if (deffile[fnm->ftp].ntps)
    {
<<<<<<< HEAD
        set_grpfnm(fnm,name ? buf : NULL,bCanNotOverride);
    }
    else
    {
        if ((name == NULL) || !(bCanNotOverride || (default_file_name ==NULL)))
        {
            const char *defnm=ftp2defnm(fnm->ftp);
            strcpy(buf,defnm);
        }
        set_extension(buf,fnm->ftp);
        
=======
        set_grpfnm(fnm, name ? buf : NULL, bCanNotOverride);
    }
    else
    {
        if ((name == NULL) || !(bCanNotOverride || (default_file_name == NULL)))
        {
            const char *defnm = ftp2defnm(fnm->ftp);
            strcpy(buf, defnm);
        }
        set_extension(buf, fnm->ftp);

>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        add_filenm(fnm, buf);
    }
}

static void set_filenms(int nf, t_filenm fnm[], gmx_bool bReadNode)
{
    int i;

    for (i = 0; (i < nf); i++)
<<<<<<< HEAD
        if (!IS_SET(fnm[i]))
            set_filenm(&(fnm[i]), fnm[i].fn, FALSE, bReadNode);
=======
    {
        if (!IS_SET(fnm[i]))
        {
            set_filenm(&(fnm[i]), fnm[i].fn, FALSE, bReadNode);
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

void parse_file_args(int *argc, char *argv[], int nf, t_filenm fnm[],
                     gmx_bool bKeep, gmx_bool bReadNode)
{
<<<<<<< HEAD
    int i, j;
=======
    int       i, j;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    gmx_bool *bRemove;

    check_opts(nf, fnm);

    for (i = 0; (i < nf); i++)
<<<<<<< HEAD
        UN_SET(fnm[i]);

    if (*argc > 1)
    {
        snew(bRemove,(*argc)+1);
=======
    {
        UN_SET(fnm[i]);
    }

    if (*argc > 1)
    {
        snew(bRemove, (*argc)+1);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        i = 1;
        do
        {
            for (j = 0; (j < nf); j++)
            {
                if (strcmp(argv[i], fnm[j].opt) == 0)
                {
                    DO_SET(fnm[j]);
                    bRemove[i] = TRUE;
                    i++;
                    /* check if we are out of arguments for this option */
                    if ((i >= *argc) || (argv[i][0] == '-'))
<<<<<<< HEAD
                        set_filenm(&fnm[j], fnm[j].fn, FALSE, bReadNode);
=======
                    {
                        set_filenm(&fnm[j], fnm[j].fn, FALSE, bReadNode);
                    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    /* sweep up all file arguments for this option */
                    while ((i < *argc) && (argv[i][0] != '-'))
                    {
                        set_filenm(&fnm[j], argv[i], TRUE, bReadNode);
                        bRemove[i] = TRUE;
                        i++;
                        /* only repeat for 'multiple' file options: */
                        if (!IS_MULT(fnm[j]))
<<<<<<< HEAD
                            break;
=======
                        {
                            break;
                        }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                    }

                    break; /* jump out of 'j' loop */
                }
            }
            /* No file found corresponding to option argv[i] */
            if (j == nf)
<<<<<<< HEAD
                i++;
        } while (i < *argc);
=======
            {
                i++;
            }
        }
        while (i < *argc);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        if (!bKeep)
        {
            /* Remove used entries */
            for (i = j = 0; (i <= *argc); i++)
            {
                if (!bRemove[i])
<<<<<<< HEAD
                    argv[j++] = argv[i];
=======
                {
                    argv[j++] = argv[i];
                }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            }
            (*argc) = j - 1;
        }
        sfree(bRemove);
    }

    set_filenms(nf, fnm, bReadNode);

}

const char *opt2fn(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
=======
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            return fnm[i].fns[0];
        }
<<<<<<< HEAD
=======
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fprintf(stderr, "No option %s\n", opt);

    return NULL;
}

const char *opt2fn_master(const char *opt, int nfile, const t_filenm fnm[],
                          t_commrec *cr)
{
    return SIMMASTER(cr) ? opt2fn(opt, nfile, fnm) : NULL;
}

int opt2fns(char **fns[], const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
=======
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            *fns = fnm[i].fns;
            return fnm[i].nfiles;
        }
<<<<<<< HEAD
=======
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fprintf(stderr, "No option %s\n", opt);
    return 0;
}

const char *ftp2fn(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
        if (ftp == fnm[i].ftp)
            return fnm[i].fns[0];
=======
    {
        if (ftp == fnm[i].ftp)
        {
            return fnm[i].fns[0];
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fprintf(stderr, "ftp2fn: No filetype %s\n", deffile[ftp].ext);
    return NULL;
}

int ftp2fns(char **fns[], int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
=======
    {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        if (ftp == fnm[i].ftp)
        {
            *fns = fnm[i].fns;
            return fnm[i].nfiles;
        }
<<<<<<< HEAD
=======
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fprintf(stderr, "ftp2fn: No filetype %s\n", deffile[ftp].ext);
    return 0;
}

gmx_bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
        if (ftp == fnm[i].ftp)
            return (gmx_bool) IS_SET(fnm[i]);
=======
    {
        if (ftp == fnm[i].ftp)
        {
            return (gmx_bool) IS_SET(fnm[i]);
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fprintf(stderr, "ftp2bSet: No filetype %s\n", deffile[ftp].ext);

    return FALSE;
}

gmx_bool opt2bSet(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
        if (strcmp(opt, fnm[i].opt) == 0)
            return (gmx_bool) IS_SET(fnm[i]);
=======
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            return (gmx_bool) IS_SET(fnm[i]);
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    fprintf(stderr, "No option %s\n", opt);

    return FALSE;
}

const char *opt2fn_null(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
                return NULL;
            else
                return fnm[i].fns[0];
        }
=======
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return NULL;
            }
            else
            {
                return fnm[i].fns[0];
            }
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    fprintf(stderr, "No option %s\n", opt);
    return NULL;
}

const char *ftp2fn_null(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
<<<<<<< HEAD
        if (ftp == fnm[i].ftp)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
                return NULL;
            else
                return fnm[i].fns[0];
        }
=======
    {
        if (ftp == fnm[i].ftp)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return NULL;
            }
            else
            {
                return fnm[i].fns[0];
            }
        }
    }
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    fprintf(stderr, "ftp2fn: No filetype %s\n", deffile[ftp].ext);
    return NULL;
}

#if 0
<<<<<<< HEAD
static void add_filters(char *filter,int *n,int nf,const int ftp[])
{
    char buf[8];
    int i;

    sprintf(filter,"*.{");
    for(i=0; (i<nf); i++)
    {
        sprintf(buf,"%s",ftp2ext(ftp[i]));
        if (*n > 0)
            strcat(filter,",");
        strcat(filter,buf);
        (*n) ++;
    }
    strcat(filter,"}");
=======
static void add_filters(char *filter, int *n, int nf, const int ftp[])
{
    char buf[8];
    int  i;

    sprintf(filter, "*.{");
    for (i = 0; (i < nf); i++)
    {
        sprintf(buf, "%s", ftp2ext(ftp[i]));
        if (*n > 0)
        {
            strcat(filter, ",");
        }
        strcat(filter, buf);
        (*n)++;
    }
    strcat(filter, "}");
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
}

char *ftp2filter(int ftp)
{
<<<<<<< HEAD
    int n;
    static char filter[128];

    filter[0] = '\0';
    n = 0;
    switch (ftp)
    {
    case efTRX:
        add_filters(filter,&n,NTRXS,trxs);
        break;
    case efTRN:
        add_filters(filter,&n,NTRNS,trns);
        break;
    case efSTO:
        add_filters(filter,&n,NSTOS,stos);
        break;
    case efSTX:
        add_filters(filter,&n,NSTXS,stxs);
        break;
    case efTPX:
        add_filters(filter,&n,NTPXS,tpxs);
        break;
    default:
        sprintf(filter,"*%s",ftp2ext(ftp));
        break;
=======
    int         n;
    static char filter[128];

    filter[0] = '\0';
    n         = 0;
    switch (ftp)
    {
        case efTRX:
            add_filters(filter, &n, NTRXS, trxs);
            break;
        case efTRN:
            add_filters(filter, &n, NTRNS, trns);
            break;
        case efSTO:
            add_filters(filter, &n, NSTOS, stos);
            break;
        case efSTX:
            add_filters(filter, &n, NSTXS, stxs);
            break;
        case efTPX:
            add_filters(filter, &n, NTPXS, tpxs);
            break;
        default:
            sprintf(filter, "*%s", ftp2ext(ftp));
            break;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    }
    return filter;
}
#endif

gmx_bool is_optional(const t_filenm *fnm)
{
    return ((fnm->flag & ffOPT) == ffOPT);
}

gmx_bool is_output(const t_filenm *fnm)
{
    return ((fnm->flag & ffWRITE) == ffWRITE);
}

gmx_bool is_set(const t_filenm *fnm)
{
    return ((fnm->flag & ffSET) == ffSET);
}

int add_suffix_to_output_names(t_filenm *fnm, int nfile, const char *suffix)
{
<<<<<<< HEAD
    int i, j, pos;
    char buf[STRLEN], newname[STRLEN];
=======
    int   i, j, pos;
    char  buf[STRLEN], newname[STRLEN];
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    char *extpos;

    for (i = 0; i < nfile; i++)
    {
        if (is_output(&fnm[i]) && fnm[i].ftp != efCPT)
        {
<<<<<<< HEAD
            /* We never use multiple _outputs_, but we might as well check 
             for it, just in case... */
            for (j = 0; j < fnm[i].nfiles; j++)
            {
                strncpy(buf, fnm[i].fns[j], STRLEN - 1);
                extpos = strrchr(buf, '.');
=======
            /* We never use multiple _outputs_, but we might as well check
               for it, just in case... */
            for (j = 0; j < fnm[i].nfiles; j++)
            {
                strncpy(buf, fnm[i].fns[j], STRLEN - 1);
                extpos  = strrchr(buf, '.');
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                *extpos = '\0';
                sprintf(newname, "%s%s.%s", buf, suffix, extpos + 1);
                free(fnm[i].fns[j]);
                fnm[i].fns[j] = strdup(newname);
            }
        }
    }
    return 0;
}

t_filenm *dup_tfn(int nf, const t_filenm tfn[])
{
<<<<<<< HEAD
    int i, j;
=======
    int       i, j;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    t_filenm *ret;

    snew(ret, nf);
    for (i = 0; i < nf; i++)
    {
        ret[i] = tfn[i]; /* just directly copy all non-string fields */
        if (tfn[i].opt)
<<<<<<< HEAD
            ret[i].opt = strdup(tfn[i].opt);
        else
            ret[i].opt = NULL;

        if (tfn[i].fn)
            ret[i].fn = strdup(tfn[i].fn);
        else
            ret[i].fn = NULL;

        if (tfn[i].nfiles > 0)
        {
            snew(ret[i].fns,tfn[i].nfiles);
=======
        {
            ret[i].opt = strdup(tfn[i].opt);
        }
        else
        {
            ret[i].opt = NULL;
        }

        if (tfn[i].fn)
        {
            ret[i].fn = strdup(tfn[i].fn);
        }
        else
        {
            ret[i].fn = NULL;
        }

        if (tfn[i].nfiles > 0)
        {
            snew(ret[i].fns, tfn[i].nfiles);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            for (j = 0; j < tfn[i].nfiles; j++)
            {
                ret[i].fns[j] = strdup(tfn[i].fns[j]);
            }
        }
    }
    return ret;
}

void done_filenms(int nf, t_filenm fnm[])
{
    int i, j;

    for (i = 0; i < nf; ++i)
    {
        for (j = 0; j < fnm[i].nfiles; ++j)
        {
            sfree(fnm[i].fns[j]);
        }
        sfree(fnm[i].fns);
        fnm[i].fns = NULL;
    }
}
