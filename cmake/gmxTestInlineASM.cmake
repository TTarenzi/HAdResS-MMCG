<<<<<<< HEAD
=======
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
#
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
# - Define macro to check GCC x86 inline ASM support
#
#  GMX_TEST_INLINE_ASM_GCC_X86(VARIABLE)
#
#  VARIABLE will be set to true if GCC x86 inline asm works.
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_INLINE_ASM_GCC_X86 VARIABLE)
    IF(NOT DEFINED ${VARIABLE})
        
        MESSAGE(STATUS "Checking for GCC x86 inline asm")

        TRY_COMPILE(${VARIABLE} "${CMAKE_BINARY_DIR}"    
<<<<<<< HEAD
                    "${CMAKE_SOURCE_DIR}/cmake/TestInlineASM_gcc_x86.c")
=======
                    "${CMAKE_SOURCE_DIR}/cmake/TestInlineASM_gcc_x86.c"
                    OUTPUT_VARIABLE INLINE_ASM_COMPILE_OUTPUT)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

        if(${VARIABLE})
            MESSAGE(STATUS "Checking for GCC x86 inline asm - supported")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for GCC x86 inline asm" FORCE)
        else(${VARIABLE})
            MESSAGE(STATUS "Checking for GCC x86 inline asm - not supported")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for GCC x86 inline asm" FORCE)
      	endif(${VARIABLE})

<<<<<<< HEAD

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_INLINE_ASM_GCC_X86 VARIABLE)



<<<<<<< HEAD
# - Define macro to check MSVC x86 inline ASM support
#
#  GMX_TEST_INLINE_ASM_MSVC_X86(VARIABLE)
#
#  VARIABLE will be set to true if MSVC x86 inline asm works.
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_INLINE_ASM_MSVC_X86 VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for MSVC x86 inline asm")

        TRY_COMPILE(${VARIABLE} "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestInlineASM_msvc_x86.c")

        if(${VARIABLE})
            MESSAGE(STATUS "Checking for MSVC x86 inline asm - supported")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for MSVC x86 inline asm" FORCE)
      	else(${VARIABLE})
            MESSAGE(STATUS "Checking for MSVC x86 inline asm - not supported")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for MSVC x86 inline asm" FORCE)
        endif(${VARIABLE})

    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_INLINE_ASM_MSVC_X86 VARIABLE)




=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
