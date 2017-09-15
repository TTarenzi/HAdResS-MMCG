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
 * Gromacs Runs On Most of All Computer Systems
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

#ifndef _writeps_h
#define _writeps_h


#include <stdio.h>
<<<<<<< HEAD
=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
typedef enum { 
  eXCenter, eXLeft, eXRight
} eXPos;

typedef enum { 
  eYCenter, eYTop,  eYBottom
} eYPos; 

enum { 
  efontTIMES, efontTIMESITALIC, efontTIMESBOLD, efontTIMESBOLDITALIC,
  efontHELV,  efontHELVITALIC,  efontHELVBOLD,  efontHELVBOLDITALIC,
  efontCOUR,  efontCOURITALIC,  efontCOURBOLD,  efontCOURBOLDITALIC,
  efontNR };
=======
typedef enum {
    eXCenter, eXLeft, eXRight
} eXPos;

typedef enum {
    eYCenter, eYTop,  eYBottom
} eYPos;

enum {
    efontTIMES, efontTIMESITALIC, efontTIMESBOLD, efontTIMESBOLDITALIC,
    efontHELV,  efontHELVITALIC,  efontHELVBOLD,  efontHELVBOLDITALIC,
    efontCOUR,  efontCOURITALIC,  efontCOURBOLD,  efontCOURBOLDITALIC,
    efontNR
};
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


typedef struct t_int_psdata *t_psdata;
/* Only use t_psdata - it is a pointer to an abstract datatype
 * that maintains the state of the postscript currently written.
 */

extern const char *fontnm[efontNR];

<<<<<<< HEAD
t_psdata ps_open(const char *fn,real x1,real y1,real x2,real y2);

void ps_linewidth(t_psdata ps,int lw);
void ps_color(t_psdata ps,real r,real g,real b);
void ps_rgb(t_psdata ps,t_rgb *rgb);

void ps_rgb_box(t_psdata ps,t_rgb *rgb);
void ps_rgb_nbox(t_psdata ps,t_rgb *rgb,real n);
void ps_init_rgb_box(t_psdata ps,real xbox, real ybox);
void ps_init_rgb_nbox(t_psdata ps,real xbox, real ybox);

void ps_lineto(t_psdata ps,real x,real y);
void ps_linerel(t_psdata ps,real dx,real dy);

void ps_moveto(t_psdata ps,real x,real y);
void ps_moverel(t_psdata ps,real dx,real dy);

void ps_line(t_psdata ps,real x1,real y1,real x2,real y2);

void ps_box(t_psdata ps,real x1,real y1,real x2,real y2);
void ps_fillbox(t_psdata ps,real x1,real y1,real x2,real y2);

void ps_arc(t_psdata ps,real x1,real y1,real rad,real a0,real a1);
void ps_fillarc(t_psdata ps,real x1,real y1,real rad,real a0,real a1);
void ps_arcslice(t_psdata ps,real xc,real yc,
			real rad1,real rad2,real a0,real a1);
void ps_fillarcslice(t_psdata ps,real xc,real yc,
			    real rad1,real rad2,real a0,real a1);

void ps_circle(t_psdata ps,real x1,real y1,real rad);

void ps_font(t_psdata ps,int font,real size);
void ps_strfont(t_psdata ps,char *font,real size);

void ps_text(t_psdata ps,real x1,real y1,const char *str);
void ps_ctext(t_psdata ps,real x1,real y1,const char *str,int expos);

void ps_close(t_psdata ps);

void ps_flip(t_psdata ps,gmx_bool bPlus);
/* Rotate over 90 (bPlus) or -90 (!bPlus) degrees */

void ps_rotate(t_psdata ps,real angle);

void ps_translate(t_psdata ps,real x,real y);

void ps_setorigin(t_psdata ps);
=======
GMX_LIBGMX_EXPORT
t_psdata ps_open(const char *fn, real x1, real y1, real x2, real y2);

GMX_LIBGMX_EXPORT
void ps_linewidth(t_psdata ps, int lw);
GMX_LIBGMX_EXPORT
void ps_color(t_psdata ps, real r, real g, real b);
GMX_LIBGMX_EXPORT
void ps_rgb(t_psdata ps, t_rgb *rgb);

void ps_rgb_box(t_psdata ps, t_rgb *rgb);
GMX_LIBGMX_EXPORT
void ps_rgb_nbox(t_psdata ps, t_rgb *rgb, real n);
GMX_LIBGMX_EXPORT
void ps_init_rgb_box(t_psdata ps, real xbox, real ybox);
GMX_LIBGMX_EXPORT
void ps_init_rgb_nbox(t_psdata ps, real xbox, real ybox);

void ps_lineto(t_psdata ps, real x, real y);
void ps_linerel(t_psdata ps, real dx, real dy);

GMX_LIBGMX_EXPORT
void ps_moveto(t_psdata ps, real x, real y);
GMX_LIBGMX_EXPORT
void ps_moverel(t_psdata ps, real dx, real dy);

GMX_LIBGMX_EXPORT
void ps_line(t_psdata ps, real x1, real y1, real x2, real y2);

GMX_LIBGMX_EXPORT
void ps_box(t_psdata ps, real x1, real y1, real x2, real y2);
GMX_LIBGMX_EXPORT
void ps_fillbox(t_psdata ps, real x1, real y1, real x2, real y2);

void ps_arc(t_psdata ps, real x1, real y1, real rad, real a0, real a1);
void ps_fillarc(t_psdata ps, real x1, real y1, real rad, real a0, real a1);
GMX_LIBGMX_EXPORT
void ps_arcslice(t_psdata ps, real xc, real yc,
                 real rad1, real rad2, real a0, real a1);
GMX_LIBGMX_EXPORT
void ps_fillarcslice(t_psdata ps, real xc, real yc,
                     real rad1, real rad2, real a0, real a1);

void ps_circle(t_psdata ps, real x1, real y1, real rad);

GMX_LIBGMX_EXPORT
void ps_font(t_psdata ps, int font, real size);
GMX_LIBGMX_EXPORT
void ps_strfont(t_psdata ps, char *font, real size);

void ps_text(t_psdata ps, real x1, real y1, const char *str);
GMX_LIBGMX_EXPORT
void ps_ctext(t_psdata ps, real x1, real y1, const char *str, int expos);

GMX_LIBGMX_EXPORT
void ps_close(t_psdata ps);

GMX_LIBGMX_EXPORT
void ps_flip(t_psdata ps, gmx_bool bPlus);
/* Rotate over 90 (bPlus) or -90 (!bPlus) degrees */

GMX_LIBGMX_EXPORT
void ps_rotate(t_psdata ps, real angle);

GMX_LIBGMX_EXPORT
void ps_translate(t_psdata ps, real x, real y);

GMX_LIBGMX_EXPORT
void ps_setorigin(t_psdata ps);
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
void ps_unsetorigin(t_psdata ps);

void viewps(char *fn);

<<<<<<< HEAD
void ps_comment(t_psdata ps,const char *s);
=======
GMX_LIBGMX_EXPORT
void ps_comment(t_psdata ps, const char *s);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

#ifdef __cplusplus
}
#endif

<<<<<<< HEAD
#endif	/* _writeps_h */

=======
#endif  /* _writeps_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
