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
 * Gyas ROwers Mature At Cryogenic Speed
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

#ifndef _pulldown_h
#define _pulldown_h

#include "popup.h"

typedef struct {
<<<<<<< HEAD
  t_windata wd;
  int       nmenu;
  int       nsel;
  int       *xpos;
  t_menu    **m;
  const char **title;
} t_pulldown;

extern t_pulldown *init_pd(t_x11 *x11,Window Parent,int width,int height,
			   unsigned long fg,unsigned long bg,
			   int nmenu,int *nsub,t_mentry *ent[],
=======
    t_windata    wd;
    int          nmenu;
    int          nsel;
    int         *xpos;
    t_menu     **m;
    const char **title;
} t_pulldown;

extern t_pulldown *init_pd(t_x11 *x11, Window Parent, int width, int height,
                           unsigned long fg, unsigned long bg,
                           int nmenu, int *nsub, t_mentry *ent[],
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                           const char **title);
/* nmenu is the number of submenus, title are the titles of
 * the submenus, nsub are the numbers of entries in each submenu
 * ent are the entries in the pulldown menu, analogous to these in the
 * popup menu.
 * The Parent is the parent window, the width is the width of the parent
 * window. The Menu is constructed as a bar on the topside of the window
 * (as usual). It calculates it own height by the font size.
 * !!!
 * !!! Do not destroy the ent structure, or the titles, while using
 * !!! the menu.
 * !!!
 * When the menu is selected, a ClientMessage will be sent to the Parent
 * specifying the selected item in xclient.data.l[0].
 */

<<<<<<< HEAD
extern void hide_pd(t_x11 *x11,t_pulldown *pd);
/* Hides any menu that is still on the screen when it shouldn't */

extern void check_pd_item(t_pulldown *pd,int nreturn,gmx_bool bStatus);
=======
extern void hide_pd(t_x11 *x11, t_pulldown *pd);
/* Hides any menu that is still on the screen when it shouldn't */

extern void check_pd_item(t_pulldown *pd, int nreturn, gmx_bool bStatus);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Set the bChecked field in the pd item with return code
 * nreturn to bStatus. This function must always be called when
 * the bChecked flag has to changed.
 */

<<<<<<< HEAD
extern void done_pd(t_x11 *x11,t_pulldown *pd);
=======
extern void done_pd(t_x11 *x11, t_pulldown *pd);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* This routine destroys the menu pd, and unregisters it with x11 */

extern int pd_width(t_pulldown *pd);
/* Return the width of the window */

extern int pd_height(t_pulldown *pd);
/* Return the height of the window */

<<<<<<< HEAD
#endif	/* _pulldown_h */
=======
#endif  /* _pulldown_h */
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2