/*
<<<<<<< HEAD
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
=======
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
/*! \internal \file
 * \brief Internal header file used by the selection tokenizer.
 */
#ifndef SELECTION_SCANNER_INTERNAL_H
#define SELECTION_SCANNER_INTERNAL_H

#include "parser.h"

/** The scanning function generated by Flex. */
#define YY_DECL int _gmx_sel_yylex(YYSTYPE *yylval, yyscan_t yyscanner)
YY_DECL;

/* These need to be defined before including scanner_flex.h, because it
 * uses YY_EXTRA_TYPE. But we also need to include it before defining
 * gmx_sel_lexer_t; hence the forward declaration. */
struct gmx_sel_lexer_t;
#define YY_EXTRA_TYPE struct gmx_sel_lexer_t *

/* We cannot include scanner_flex.h from the scanner itself, because it
 * seems to break everything. */
/* And we need to define YY_NO_UNISTD_H here as well, otherwise unistd.h
 * gets included in other files than scanner.c... */
#ifndef FLEX_SCANNER
#define YY_NO_UNISTD_H
#include "scanner_flex.h"
#endif

/*! \brief
 * Internal data structure for the selection tokenizer state.
 */
typedef struct gmx_sel_lexer_t
{
<<<<<<< HEAD
    struct gmx_ana_selcollection_t  *sc;
    struct gmx_ana_indexgrps_t      *grps;
    int                              nexpsel;

    gmx_bool                             bInteractive;
    char                            *inputstr;
    int                              nalloc_input;

    char                            *pselstr;
    int                              pslen;
    int                              nalloc_psel;

    struct gmx_ana_selmethod_t     **mstack;
    int                              msp;
    int                              mstack_alloc;

    int                              neom;
    struct gmx_ana_selparam_t       *nextparam;
    gmx_bool                             bBoolNo;
    struct gmx_ana_selmethod_t      *nextmethod;
    int                              prev_pos_kw;
=======
    struct gmx_ana_selcollection_t      *sc;
    struct gmx_ana_indexgrps_t          *grps;
    int                                  nexpsel;

    gmx_bool                             bInteractive;
    char                                *inputstr;
    int                                  nalloc_input;

    char                                *pselstr;
    int                                  pslen;
    int                                  nalloc_psel;

    struct gmx_ana_selmethod_t         **mstack;
    int                                  msp;
    int                                  mstack_alloc;

    int                                  neom;
    struct gmx_ana_selparam_t           *nextparam;
    gmx_bool                             bBoolNo;
    struct gmx_ana_selmethod_t          *nextmethod;
    int                                  prev_pos_kw;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

    gmx_bool                             bMatchOf;
    gmx_bool                             bMatchBool;
    gmx_bool                             bCmdStart;

    gmx_bool                             bBuffer;
<<<<<<< HEAD
    YY_BUFFER_STATE                  buffer;
=======
    YY_BUFFER_STATE                      buffer;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
} gmx_sel_lexer_t;

/* Because Flex defines yylval, yytext, and yyleng as macros,
 * and this file is included from scanner.l,
 * we cannot have them here as parameter names... */
/** Internal function for cases where several tokens need to be returned. */
int
_gmx_sel_lexer_process_pending(YYSTYPE *, gmx_sel_lexer_t *state);
/** Internal function that processes identifier tokens. */
int
<<<<<<< HEAD
_gmx_sel_lexer_process_identifier(YYSTYPE *, char *, size_t,
                                  gmx_sel_lexer_t *state);
=======
    _gmx_sel_lexer_process_identifier(YYSTYPE *, char *, size_t,
                                      gmx_sel_lexer_t *state);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
/** Internal function to add a token to the pretty-printed selection text. */
void
_gmx_sel_lexer_add_token(const char *str, int len, gmx_sel_lexer_t *state);

#endif
