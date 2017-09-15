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

#ifndef _3dview_h
#define _3dview_h
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

#define WW 3

typedef real vec4[4];

typedef real mat4[4][4];

typedef int  iv2[2];
<<<<<<< HEAD
 
typedef struct {
  matrix box;
  int    ecenter;       /* enum for centering, see pbc.h        */
  vec4   eye,origin;	/* The eye and origin position		*/
  mat4   proj;		/* Projection matrix 			*/
  mat4   Rot;           /* Total rotation matrix                */
  real   sc_x,sc_y;	/* Scaling for aspect ratio		*/
  mat4   RotP[DIM];     /* state for 3d rotations  */
  mat4   RotM[DIM];
} t_3dview;

extern void print_m4(FILE *fp,const char *s,mat4 A);

extern void print_v4(FILE *fp,char *s,int dim,real *a);

extern void m4_op(mat4 m,rvec x,vec4 v);

extern void unity_m4(mat4 m);

extern void mult_matrix(mat4 A, mat4 B, mat4 C);

extern void rotate(int axis, real angle, mat4 A);

extern void translate(real tx, real ty, real tz, mat4 A);

extern void m4_op(mat4 m,rvec x,vec4 v);
=======

typedef struct {
    matrix box;
    int    ecenter;     /* enum for centering, see pbc.h        */
    vec4   eye, origin; /* The eye and origin position		*/
    mat4   proj;        /* Projection matrix            */
    mat4   Rot;         /* Total rotation matrix                */
    real   sc_x, sc_y;  /* Scaling for aspect ratio		*/
    mat4   RotP[DIM];   /* state for 3d rotations  */
    mat4   RotM[DIM];
} t_3dview;

GMX_LIBGMX_EXPORT
extern void print_m4(FILE *fp, const char *s, mat4 A);

extern void print_v4(FILE *fp, char *s, int dim, real *a);

GMX_LIBGMX_EXPORT
extern void m4_op(mat4 m, rvec x, vec4 v);

extern void unity_m4(mat4 m);

GMX_LIBGMX_EXPORT
extern void mult_matrix(mat4 A, mat4 B, mat4 C);

GMX_LIBGMX_EXPORT
extern void rotate(int axis, real angle, mat4 A);

GMX_LIBGMX_EXPORT
extern void translate(real tx, real ty, real tz, mat4 A);

GMX_LIBGMX_EXPORT
extern void m4_op(mat4 m, rvec x, vec4 v);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

extern void calculate_view(t_3dview *view);

extern t_3dview *init_view(matrix box);
/* Generate the view matrix from the eye pos and the origin,
 * applying also the scaling for the aspect ration.
 * There is no accompanying done_view routine: the struct can simply
 * be sfree'd.
 */

/* The following options are present on the 3d struct:
 * zoom (scaling)
 * rotate around the center of the box
 * reset the view
 */

<<<<<<< HEAD
extern gmx_bool zoom_3d(t_3dview *view,real fac);
=======
extern gmx_bool zoom_3d(t_3dview *view, real fac);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Zoom in or out with factor fac, returns TRUE when zoom successful,
 * FALSE otherwise.
 */

extern void init_rotate_3d(t_3dview *view);
/* Initiates the state of 3d rotation matrices in the structure */

<<<<<<< HEAD
extern void rotate_3d(t_3dview *view,int axis,gmx_bool bPositive);
/* Rotate the eye around the center of the box, around axis */

extern void translate_view(t_3dview *view,int axis,gmx_bool bPositive);
=======
extern void rotate_3d(t_3dview *view, int axis, gmx_bool bPositive);
/* Rotate the eye around the center of the box, around axis */

extern void translate_view(t_3dview *view, int axis, gmx_bool bPositive);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/* Translate the origin at which one is looking */

extern void reset_view(t_3dview *view);
/* Reset the viewing to the initial view */

#ifdef __cplusplus
}
#endif

#endif
<<<<<<< HEAD

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2