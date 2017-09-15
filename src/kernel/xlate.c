<<<<<<< HEAD
/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
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
#include <string.h>
#include "typedefs.h"
#include "strdb.h"
#include "string2.h"
#include "smalloc.h"
#include "symtab.h"
#include "index.h"
#include "futil.h"
#include "fflibutil.h"
#include "hackblock.h"
#include "gmx_fatal.h"
#include "xlate.h"

typedef struct {
    char *filebase;
    char *res;
    char *atom;
    char *replace;
} t_xlate_atom;

<<<<<<< HEAD
static void get_xlatoms(const char *fn,FILE *fp,
                        int *nptr,t_xlate_atom **xlptr)
{
    char filebase[STRLEN];
    char line[STRLEN];
    char abuf[1024],rbuf[1024],repbuf[1024],dumbuf[1024];
    char *_ptr;
    int  n,na,idum;
    t_xlate_atom *xl;

    fflib_filename_base(fn,filebase,STRLEN);
=======
static void get_xlatoms(const char *fn, FILE *fp,
                        int *nptr, t_xlate_atom **xlptr)
{
    char          filebase[STRLEN];
    char          line[STRLEN];
    char          abuf[1024], rbuf[1024], repbuf[1024], dumbuf[1024];
    char         *_ptr;
    int           n, na, idum;
    t_xlate_atom *xl;

    fflib_filename_base(fn, filebase, STRLEN);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    n  = *nptr;
    xl = *xlptr;

<<<<<<< HEAD
    while (get_a_line(fp,line,STRLEN))
    {
        na = sscanf(line,"%s%s%s%s",rbuf,abuf,repbuf,dumbuf);
        /* Check if we are reading an old format file with the number of items
         * on the first line.
         */
        if (na == 1 && n == *nptr && sscanf(rbuf,"%d",&idum) == 1)
=======
    while (get_a_line(fp, line, STRLEN))
    {
        na = sscanf(line, "%s%s%s%s", rbuf, abuf, repbuf, dumbuf);
        /* Check if we are reading an old format file with the number of items
         * on the first line.
         */
        if (na == 1 && n == *nptr && sscanf(rbuf, "%d", &idum) == 1)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            continue;
        }
        if (na != 3)
        {
<<<<<<< HEAD
            gmx_fatal(FARGS,"Expected a residue name and two atom names in file '%s', not '%s'",fn,line);
        }
        
        srenew(xl,n+1);
        xl[n].filebase = strdup(filebase);

        /* Use wildcards... */
        if (strcmp(rbuf,"*") != 0)
=======
            gmx_fatal(FARGS, "Expected a residue name and two atom names in file '%s', not '%s'", fn, line);
        }

        srenew(xl, n+1);
        xl[n].filebase = strdup(filebase);

        /* Use wildcards... */
        if (strcmp(rbuf, "*") != 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        {
            xl[n].res = strdup(rbuf);
        }
        else
        {
            xl[n].res = NULL;
        }
<<<<<<< HEAD
        
        /* Replace underscores in the string by spaces */
        while ((_ptr = strchr(abuf,'_')) != 0)
        {
            *_ptr = ' ';
        }
        
        xl[n].atom = strdup(abuf);
=======

        /* Replace underscores in the string by spaces */
        while ((_ptr = strchr(abuf, '_')) != 0)
        {
            *_ptr = ' ';
        }

        xl[n].atom    = strdup(abuf);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        xl[n].replace = strdup(repbuf);
        n++;
    }

    *nptr  = n;
    *xlptr = xl;
}

<<<<<<< HEAD
static void done_xlatom(int nxlate,t_xlate_atom *xlatom)
{
    int i;
    
    for(i=0; (i<nxlate); i++)
=======
static void done_xlatom(int nxlate, t_xlate_atom *xlatom)
{
    int i;

    for (i = 0; (i < nxlate); i++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        sfree(xlatom[i].filebase);
        if (xlatom[i].res != NULL)
        {
            sfree(xlatom[i].res);
        }
        sfree(xlatom[i].atom);
        sfree(xlatom[i].replace);
    }
    sfree(xlatom);
}

<<<<<<< HEAD
void rename_atoms(const char *xlfile,const char *ffdir,
                  t_atoms *atoms,t_symtab *symtab,const t_restp *restp,
                  gmx_bool bResname,gmx_residuetype_t rt,gmx_bool bReorderNum,
                  gmx_bool bVerbose)
{
    FILE *fp;
    int nxlate,a,i,resind;
    t_xlate_atom *xlatom;
    int  nf;
    char **f;
    char c,*rnm,atombuf[32],*ptr0,*ptr1;
    gmx_bool bReorderedNum,bRenamed,bMatch;
=======
void rename_atoms(const char *xlfile, const char *ffdir,
                  t_atoms *atoms, t_symtab *symtab, const t_restp *restp,
                  gmx_bool bResname, gmx_residuetype_t rt, gmx_bool bReorderNum,
                  gmx_bool bVerbose)
{
    FILE         *fp;
    int           nxlate, a, i, resind;
    t_xlate_atom *xlatom;
    int           nf;
    char        **f;
    char          c, *rnm, atombuf[32], *ptr0, *ptr1;
    gmx_bool      bReorderedNum, bRenamed, bMatch;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    nxlate = 0;
    xlatom = NULL;
    if (xlfile != NULL)
    {
        fp = libopen(xlfile);
<<<<<<< HEAD
        get_xlatoms(xlfile,fp,&nxlate,&xlatom);
=======
        get_xlatoms(xlfile, fp, &nxlate, &xlatom);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        fclose(fp);
    }
    else
    {
<<<<<<< HEAD
        nf = fflib_search_file_end(ffdir,".arn",FALSE,&f);
        for(i=0; i<nf; i++)
        {
            fp = fflib_open(f[i]);
            get_xlatoms(f[i],fp,&nxlate,&xlatom);
            fclose(fp);
=======
        nf = fflib_search_file_end(ffdir, ".arn", FALSE, &f);
        for (i = 0; i < nf; i++)
        {
            fp = fflib_open(f[i]);
            get_xlatoms(f[i], fp, &nxlate, &xlatom);
            ffclose(fp);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
            sfree(f[i]);
        }
        sfree(f);
    }

<<<<<<< HEAD
    for(a=0; (a<atoms->nr); a++)
=======
    for (a = 0; (a < atoms->nr); a++)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    {
        resind = atoms->atom[a].resind;
        if (bResname)
        {
            rnm = *(atoms->resinfo[resind].name);
        }
        else
        {
            rnm = *(atoms->resinfo[resind].rtp);
        }
<<<<<<< HEAD
               
        strcpy(atombuf,*(atoms->atomname[a]));
=======

        strcpy(atombuf, *(atoms->atomname[a]));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        bReorderedNum = FALSE;
        if (bReorderNum)
        {
            if (isdigit(atombuf[0]))
            {
                c = atombuf[0];
<<<<<<< HEAD
                for (i=0; ((size_t)i<strlen(atombuf)-1); i++)
                {
                    atombuf[i] = atombuf[i+1];
                }
                atombuf[i] = c;
=======
                for (i = 0; ((size_t)i < strlen(atombuf)-1); i++)
                {
                    atombuf[i] = atombuf[i+1];
                }
                atombuf[i]    = c;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                bReorderedNum = TRUE;
            }
        }
        bRenamed = FALSE;
<<<<<<< HEAD
        for(i=0; (i<nxlate) && !bRenamed; i++) {
            /* Check if the base file name of the rtp and arn entry match */
            if (restp == NULL ||
                gmx_strcasecmp(restp[resind].filebase,xlatom[i].filebase) == 0)
            {
                /* Match the residue name */
                bMatch = (xlatom[i].res == NULL ||
                          (gmx_strcasecmp("protein",xlatom[i].res) == 0 &&
                           gmx_residuetype_is_protein(rt,rnm)) ||
                          (gmx_strcasecmp("DNA",xlatom[i].res) == 0 &&
                           gmx_residuetype_is_dna(rt,rnm)) ||
                          (gmx_strcasecmp("RNA",xlatom[i].res) == 0 &&
                           gmx_residuetype_is_rna(rt,rnm)));
=======
        for (i = 0; (i < nxlate) && !bRenamed; i++)
        {
            /* Check if the base file name of the rtp and arn entry match */
            if (restp == NULL ||
                gmx_strcasecmp(restp[resind].filebase, xlatom[i].filebase) == 0)
            {
                /* Match the residue name */
                bMatch = (xlatom[i].res == NULL ||
                          (gmx_strcasecmp("protein", xlatom[i].res) == 0 &&
                           gmx_residuetype_is_protein(rt, rnm)) ||
                          (gmx_strcasecmp("DNA", xlatom[i].res) == 0 &&
                           gmx_residuetype_is_dna(rt, rnm)) ||
                          (gmx_strcasecmp("RNA", xlatom[i].res) == 0 &&
                           gmx_residuetype_is_rna(rt, rnm)));
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                if (!bMatch)
                {
                    ptr0 = rnm;
                    ptr1 = xlatom[i].res;
                    while (ptr0[0] != '\0' && ptr1[0] != '\0' &&
                           (ptr0[0] == ptr1[0] || ptr1[0] == '?'))
                    {
                        ptr0++;
                        ptr1++;
                    }
                    bMatch = (ptr0[0] == '\0' && ptr1[0] == '\0');
                }
<<<<<<< HEAD
                if (bMatch && strcmp(atombuf,xlatom[i].atom) == 0)
=======
                if (bMatch && strcmp(atombuf, xlatom[i].atom) == 0)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                {
                    /* We have a match. */
                    /* Don't free the old atomname,
                     * since it might be in the symtab.
                     */
                    ptr0 = strdup(xlatom[i].replace);
                    if (bVerbose)
                    {
                        printf("Renaming atom '%s' in residue %d %s to '%s'\n",
                               *atoms->atomname[a],
                               atoms->resinfo[resind].nr,
                               *atoms->resinfo[resind].name,
                               ptr0);
                    }
<<<<<<< HEAD
                    atoms->atomname[a] = put_symtab(symtab,ptr0);
                    bRenamed = TRUE;
=======
                    atoms->atomname[a] = put_symtab(symtab, ptr0);
                    bRenamed           = TRUE;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                }
            }
        }
        if (bReorderedNum && !bRenamed)
        {
<<<<<<< HEAD
            atoms->atomname[a] = put_symtab(symtab,atombuf);
        }
    }

    done_xlatom(nxlate,xlatom);
}

=======
            atoms->atomname[a] = put_symtab(symtab, atombuf);
        }
    }

    done_xlatom(nxlate, xlatom);
}
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
