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
 * \brief Implementation of functions in selhelp.c.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <macros.h>
#include <string2.h>
#include <wman.h>

#include "selcollection.h"
#include "selmethod.h"
#include "selhelp.h"
#include "symrec.h"

typedef struct {
    const char  *topic;
    int          nl;
    const char **text;
} t_selection_help_item;

<<<<<<< HEAD
static const char *help_common[] = {
=======
static const char                 *help_common[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION HELP[PAR]",

    "This program supports selections in addition to traditional index files.",
    "Please read the subtopic pages (available through \"help topic\") for",
    "more information.",
    "Explanation of command-line arguments for specifying selections can be",
    "found under the \"cmdline\" subtopic, and general selection syntax is",
    "described under \"syntax\". Available keywords can be found under",
    "\"keywords\", and concrete examples under \"examples\".",
    "Other subtopics give more details on certain aspects.",
    "\"help all\" prints the help for all subtopics.",
};

<<<<<<< HEAD
static const char *help_arithmetic[] = {
=======
static const char                 *help_arithmetic[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "ARITHMETIC EXPRESSIONS IN SELECTIONS[PAR]",

    "Basic arithmetic evaluation is supported for numeric expressions.",
    "Supported operations are addition, subtraction, negation, multiplication,",
    "division, and exponentiation (using ^).",
    "Result of a division by zero or other illegal operations is undefined.",
};

<<<<<<< HEAD
static const char *help_cmdline[] = {
=======
static const char                 *help_cmdline[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION COMMAND-LINE ARGUMENTS[PAR]",

    "There are two alternative command-line arguments for specifying",
    "selections:[BR]",
    "1. [TT]-select[tt] can be used to specify the complete selection as a",
    "string on the command line.[BR]",
    "2. [TT]-sf[tt] can be used to specify a file name from which the",
    "selection is read.[BR]",
    "If both options are specified, [TT]-select[tt] takes precedence.",
    "If neither of the above is present, the user is prompted to type the",
    "selection on the standard input (a pipe can also be used to provide",
    "the selections in this case).",
    "This is also done if an empty string is passed to [TT]-select[tt].[PAR]",

    "Option [TT]-n[tt] can be used to provide an index file.",
    "If no index file is provided, default groups are generated.",
    "In both cases, the user can also select an index group instead of",
    "writing a full selection.",
    "The default groups are generated by reading selections from a file",
    "[TT]defselection.dat[tt]. If such a file is found in the current",
    "directory, it is used instead of the one provided by default.[PAR]",

    "Depending on the tool, two additional command-line arguments may be",
    "available to control the behavior:[BR]",
    "1. [TT]-seltype[tt] can be used to specify the default type of",
    "positions to calculate for each selection.[BR]",
    "2. [TT]-selrpos[tt] can be used to specify the default type of",
    "positions used in selecting atoms by coordinates.[BR]",
    "See \"help positions\" for more information on these options.",
};

<<<<<<< HEAD
static const char *help_eval[] = {
=======
static const char                 *help_eval[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION EVALUATION AND OPTIMIZATION[PAR]",

    "Boolean evaluation proceeds from left to right and is short-circuiting",
    "i.e., as soon as it is known whether an atom will be selected, the",
    "remaining expressions are not evaluated at all.",
    "This can be used to optimize the selections: you should write the",
    "most restrictive and/or the most inexpensive expressions first in",
    "boolean expressions.",
    "The relative ordering between dynamic and static expressions does not",
    "matter: all static expressions are evaluated only once, before the first",
    "frame, and the result becomes the leftmost expression.[PAR]",

    "Another point for optimization is in common subexpressions: they are not",
    "automatically recognized, but can be manually optimized by the use of",
    "variables. This can have a big impact on the performance of complex",
    "selections, in particular if you define several index groups like this:",
    "  [TT]rdist = distance from com of resnr 1 to 5;[tt][BR]",
    "  [TT]resname RES and rdist < 2;[tt][BR]",
    "  [TT]resname RES and rdist < 4;[tt][BR]",
    "  [TT]resname RES and rdist < 6;[tt][BR]",
    "Without the variable assignment, the distances would be evaluated three",
    "times, although they are exactly the same within each selection.",
    "Anything assigned into a variable becomes a common subexpression that",
    "is evaluated only once during a frame.",
    "Currently, in some cases the use of variables can actually lead to a small",
    "performance loss because of the checks necessary to determine for which",
    "atoms the expression has already been evaluated, but this should not be",
    "a major problem.",
};

<<<<<<< HEAD
static const char *help_examples[] = {
=======
static const char                 *help_examples[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION EXAMPLES[PAR]",

    "Below, examples of increasingly complex selections are given.[PAR]",

    "Selection of all water oxygens:[BR]",
    "  resname SOL and name OW",
    "[PAR]",

    "Centers of mass of residues 1 to 5 and 10:[BR]",
    "  res_com of resnr 1 to 5 10",
    "[PAR]",

    "All atoms farther than 1 nm of a fixed position:[BR]",
    "  not within 1 of (1.2, 3.1, 2.4)",
    "[PAR]",

    "All atoms of a residue LIG within 0.5 nm of a protein (with a custom name):[BR]",
    "  \"Close to protein\" resname LIG and within 0.5 of group \"Protein\"",
    "[PAR]",

    "All protein residues that have at least one atom within 0.5 nm of a residue LIG:[BR]",
    "  group \"Protein\" and same residue as within 0.5 of resname LIG",
    "[PAR]",

    "All RES residues whose COM is between 2 and 4 nm from the COM of all of them:[BR]",
    "  rdist = res_com distance from com of resname RES[BR]",
    "  resname RES and rdist >= 2 and rdist <= 4",
    "[PAR]",

    "Selection like C1 C2 C2 C3 C3 C4 ... C8 C9 (e.g., for g_bond):[BR]",
    "  name \"C[1-8]\" merge name \"C[2-9]\"",
};

<<<<<<< HEAD
static const char *help_keywords[] = {
=======
static const char                 *help_keywords[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION KEYWORDS[PAR]",

    "The following selection keywords are currently available.",
    "For keywords marked with a star, additional help is available through",
    "\"help KEYWORD\", where KEYWORD is the name of the keyword.",
};

<<<<<<< HEAD
static const char *help_limits[] = {
=======
static const char                 *help_limits[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION LIMITATIONS[PAR]",

    "Some analysis programs may require a special structure for the input",
    "selections (e.g., [TT]g_angle[tt] requires the index group to be made",
    "of groups of three or four atoms).",
    "For such programs, it is up to the user to provide a proper selection",
    "expression that always returns such positions.",
    "[PAR]",

    "Due to technical reasons, having a negative value as the first value in",
    "expressions like[BR]",
    "[TT]charge -1 to -0.7[tt][BR]",
    "result in a syntax error. A workaround is to write[BR]",
    "[TT]charge {-1 to -0.7}[tt][BR]",
<<<<<<< HEAD
    "instead.",
};

static const char *help_positions[] = {
=======
    "instead.[PAR]",

    "When [TT]name[tt] selection keyword is used together with PDB input",
    "files, the behavior may be unintuitive. When Gromacs reads in a PDB",
    "file, 4 character atom names that start with a digit are transformed",
    "such that, e.g., 1HG2 becomes HG21, and the latter is what is matched",
    "by the [TT]name[tt] keyword. Use [TT]pdbname[tt] to match the atom name",
    "as it appears in the input PDB file.",
};

static const char                 *help_positions[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SPECIFYING POSITIONS[PAR]",

    "Possible ways of specifying positions in selections are:[PAR]",

    "1. A constant position can be defined as [TT][XX, YY, ZZ][tt], where",
    "[TT]XX[tt], [TT]YY[tt] and [TT]ZZ[tt] are real numbers.[PAR]",

    "2. [TT]com of ATOM_EXPR [pbc][tt] or [TT]cog of ATOM_EXPR [pbc][tt]",
    "calculate the center of mass/geometry of [TT]ATOM_EXPR[tt]. If",
    "[TT]pbc[tt] is specified, the center is calculated iteratively to try",
    "to deal with cases where [TT]ATOM_EXPR[tt] wraps around periodic",
    "boundary conditions.[PAR]",

    "3. [TT]POSTYPE of ATOM_EXPR[tt] calculates the specified positions for",
    "the atoms in [TT]ATOM_EXPR[tt].",
    "[TT]POSTYPE[tt] can be [TT]atom[tt], [TT]res_com[tt], [TT]res_cog[tt],",
    "[TT]mol_com[tt] or [TT]mol_cog[tt], with an optional prefix [TT]whole_[tt]",
    "[TT]part_[tt] or [TT]dyn_[tt].",
    "[TT]whole_[tt] calculates the centers for the whole residue/molecule,",
    "even if only part of it is selected.",
    "[TT]part_[tt] prefix calculates the centers for the selected atoms, but",
    "uses always the same atoms for the same residue/molecule. The used atoms",
    "are determined from the the largest group allowed by the selection.",
    "[TT]dyn_[tt] calculates the centers strictly only for the selected atoms.",
    "If no prefix is specified, whole selections default to [TT]part_[tt] and",
    "other places default to [TT]whole_[tt].",
    "The latter is often desirable to select the same molecules in different",
    "tools, while the first is a compromise between speed ([TT]dyn_[tt]",
    "positions can be slower to evaluate than [TT]part_[tt]) and intuitive",
    "behavior.[PAR]",

    "4. [TT]ATOM_EXPR[tt], when given for whole selections, is handled as 3.",
    "above, using the position type from the command-line argument",
    "[TT]-seltype[tt].[PAR]",

    "Selection keywords that select atoms based on their positions, such as",
    "[TT]dist from[tt], use by default the positions defined by the",
    "[TT]-selrpos[tt] command-line option.",
    "This can be overridden by prepending a [TT]POSTYPE[tt] specifier to the",
    "keyword. For example, [TT]res_com dist from POS[tt] evaluates the",
    "residue center of mass distances. In the example, all atoms of a residue",
    "are either selected or not, based on the single distance calculated.",
};

<<<<<<< HEAD
static const char *help_syntax[] = {
=======
static const char                 *help_syntax[] = {
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    "SELECTION SYNTAX[PAR]",

    "A set of selections consists of one or more selections, separated by",
    "semicolons. Each selection defines a set of positions for the analysis.",
    "Each selection can also be preceded by a string that gives a name for",
    "the selection for use in, e.g., graph legends.",
    "If no name is provided, the string used for the selection is used",
    "automatically as the name.[PAR]",

    "For interactive input, the syntax is slightly altered: line breaks can",
    "also be used to separate selections. \\ followed by a line break can",
    "be used to continue a line if necessary.",
    "Notice that the above only applies to real interactive input,",
    "not if you provide the selections, e.g., from a pipe.[PAR]",

    "It is possible to use variables to store selection expressions.",
    "A variable is defined with the following syntax:[BR]",
    "[TT]VARNAME = EXPR ;[tt][BR]",
    "where [TT]EXPR[tt] is any valid selection expression.",
    "After this, [TT]VARNAME[tt] can be used anywhere where [TT]EXPR[tt]",
    "would be valid.[PAR]",

    "Selections are composed of three main types of expressions, those that",
    "define atoms ([TT]ATOM_EXPR[tt]s), those that define positions",
    "([TT]POS_EXPR[tt]s), and those that evaluate to numeric values",
    "([TT]NUM_EXPR[tt]s). Each selection should be a [TT]POS_EXPR[tt]",
    "or a [TT]ATOM_EXPR[tt] (the latter is automatically converted to",
    "positions). The basic rules are as follows:[BR]",
    "1. An expression like [TT]NUM_EXPR1 < NUM_EXPR2[tt] evaluates to an",
    "[TT]ATOM_EXPR[tt] that selects all the atoms for which the comparison",
    "is true.[BR]",
    "2. Atom expressions can be combined with gmx_boolean operations such as",
    "[TT]not ATOM_EXPR[tt], [TT]ATOM_EXPR and ATOM_EXPR[tt], or",
    "[TT]ATOM_EXPR or ATOM_EXPR[tt]. Parentheses can be used to alter the",
    "evaluation order.[BR]",
    "3. [TT]ATOM_EXPR[tt] expressions can be converted into [TT]POS_EXPR[tt]",
    "expressions in various ways, see \"help positions\" for more details.[PAR]",

    "Some keywords select atoms based on string values such as the atom name.",
    "For these keywords, it is possible to use wildcards ([TT]name \"C*\"[tt])",
    "or regular expressions (e.g., [TT]resname \"R[AB]\"[tt]).",
    "The match type is automatically guessed from the string: if it contains",
    "other characters than letters, numbers, '*', or '?', it is interpreted",
    "as a regular expression.",
    "Strings that contain non-alphanumeric characters should be enclosed in",
    "double quotes as in the examples. For other strings, the quotes are",
    "optional, but if the value conflicts with a reserved keyword, a syntax",
    "error will occur. If your strings contain uppercase letters, this should",
    "not happen.[PAR]",

    "Index groups provided with the [TT]-n[tt] command-line option or",
    "generated by default can be accessed with [TT]group NR[tt] or",
    "[TT]group NAME[tt], where [TT]NR[tt] is a zero-based index of the group",
    "and [TT]NAME[tt] is part of the name of the desired group.",
    "The keyword [TT]group[tt] is optional if the whole selection is",
    "provided from an index group.",
    "To see a list of available groups in the interactive mode, press enter",
    "in the beginning of a line.",
};

static const t_selection_help_item helpitems[] = {
    {NULL,          asize(help_common),     help_common},
    {"cmdline",     asize(help_cmdline),    help_cmdline},
    {"syntax",      asize(help_syntax),     help_syntax},
    {"positions",   asize(help_positions),  help_positions},
    {"arithmetic",  asize(help_arithmetic), help_arithmetic},
    {"keywords",    asize(help_keywords),   help_keywords},
    {"evaluation",  asize(help_eval),       help_eval},
    {"limitations", asize(help_limits),     help_limits},
    {"examples",    asize(help_examples),   help_examples},
};

/*! \brief
 * Prints a brief list of keywords (selection methods) available.
 *
 * \param[in] sc    Selection collection for which the list should be printed.
 * \param[in] type  Only methods that return this type are printed.
 * \param[in] bMod  If FALSE, \ref SMETH_MODIFIER methods are excluded, otherwise
 *     only them are printed.
 */
static void
print_keyword_list(struct gmx_ana_selcollection_t *sc, e_selvalue_t type,
                   gmx_bool bMod)
{
    gmx_sel_symrec_t *symbol;

    symbol = _gmx_sel_first_symbol(sc->symtab, SYMBOL_METHOD);
    while (symbol)
    {
<<<<<<< HEAD
        gmx_ana_selmethod_t *method = _gmx_sel_sym_value_method(symbol);
=======
        gmx_ana_selmethod_t     *method = _gmx_sel_sym_value_method(symbol);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
        gmx_bool                 bShow;
        bShow = (method->type == type)
            && ((bMod && (method->flags & SMETH_MODIFIER))
                || (!bMod && !(method->flags & SMETH_MODIFIER)));
        if (bShow)
        {
            fprintf(stderr, " %c ",
                    (method->help.nlhelp > 0 && method->help.help) ? '*' : ' ');
            if (method->help.syntax)
            {
                fprintf(stderr, "%s\n", method->help.syntax);
            }
            else
            {
                const char *symname = _gmx_sel_sym_name(symbol);

                fprintf(stderr, "%s", symname);
                if (strcmp(symname, method->name) != 0)
                {
                    fprintf(stderr, " (synonym for %s)", method->name);
                }
                fprintf(stderr, "\n");
            }
        }
        symbol = _gmx_sel_next_symbol(symbol, SYMBOL_METHOD);
    }
}

/*!
 * \param[in]  sc    Selection collection for which help should be printed.
 * \param[in]  topic Topic to print help on, or NULL for general help.
 *
 * \p sc is used to get information on which keywords are available in the
 * present context.
 */
void
_gmx_sel_print_help(struct gmx_ana_selcollection_t *sc, const char *topic)
{
    const t_selection_help_item *item = NULL;
    size_t i;

    /* Find the item for the topic */
    if (!topic)
    {
        item = &helpitems[0];
    }
    else if (strcmp(topic, "all") == 0)
    {
        for (i = 0; i < asize(helpitems); ++i)
        {
            item = &helpitems[i];
            _gmx_sel_print_help(sc, item->topic);
            if (i != asize(helpitems) - 1)
            {
                fprintf(stderr, "\n\n");
            }
        }
        return;
    }
    else
    {
        for (i = 1; i < asize(helpitems); ++i)
        {
            if (strncmp(helpitems[i].topic, topic, strlen(topic)) == 0)
            {
                item = &helpitems[i];
                break;
            }
        }
    }
    /* If the topic is not found, check the available methods.
     * If they don't provide any help either, tell the user and exit. */
    if (!item)
    {
        gmx_sel_symrec_t *symbol;

        symbol = _gmx_sel_first_symbol(sc->symtab, SYMBOL_METHOD);
        while (symbol)
        {
            gmx_ana_selmethod_t *method = _gmx_sel_sym_value_method(symbol);
            if (method->help.nlhelp > 0 && method->help.help
                && strncmp(method->name, topic, strlen(topic)) == 0)
            {
                print_tty_formatted(stderr, method->help.nlhelp,
<<<<<<< HEAD
                        method->help.help, 0, NULL, NULL, FALSE);
=======
                                    method->help.help, 0, NULL, NULL, FALSE);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
                return;
            }
            symbol = _gmx_sel_next_symbol(symbol, SYMBOL_METHOD);
        }

        fprintf(stderr, "No help available for '%s'.\n", topic);
        return;
    }
    /* Print the help */
    print_tty_formatted(stderr, item->nl, item->text, 0, NULL, NULL, FALSE);
    /* Special handling of certain pages */
    if (!topic)
    {
        int len = 0;

        /* Print the subtopics on the main page */
        fprintf(stderr, "\nAvailable subtopics:\n");
        for (i = 1; i < asize(helpitems); ++i)
        {
            int len1 = strlen(helpitems[i].topic) + 2;

            len += len1;
            if (len > 79)
            {
                fprintf(stderr, "\n");
                len = len1;
            }
            fprintf(stderr, "  %s", helpitems[i].topic);
        }
        fprintf(stderr, "\n");
    }
    else if (strcmp(item->topic, "keywords") == 0)
    {
        /* Print the list of keywords */
        fprintf(stderr, "\nKeywords that select atoms by an integer property:\n");
        fprintf(stderr, "(use in expressions or like \"atomnr 1 to 5 7 9\")\n");
        print_keyword_list(sc, INT_VALUE, FALSE);

        fprintf(stderr, "\nKeywords that select atoms by a numeric property:\n");
        fprintf(stderr, "(use in expressions or like \"occupancy 0.5 to 1\")\n");
        print_keyword_list(sc, REAL_VALUE, FALSE);

        fprintf(stderr, "\nKeywords that select atoms by a string property:\n");
        fprintf(stderr, "(use like \"name PATTERN [PATTERN] ...\")\n");
        print_keyword_list(sc, STR_VALUE, FALSE);

        fprintf(stderr, "\nAdditional keywords that directly select atoms:\n");
        print_keyword_list(sc, GROUP_VALUE, FALSE);

        fprintf(stderr, "\nKeywords that directly evaluate to positions:\n");
        fprintf(stderr, "(see also \"help positions\")\n");
        print_keyword_list(sc, POS_VALUE, FALSE);

        fprintf(stderr, "\nAdditional keywords:\n");
        print_keyword_list(sc, POS_VALUE, TRUE);
        print_keyword_list(sc, NO_VALUE, TRUE);
    }
}
