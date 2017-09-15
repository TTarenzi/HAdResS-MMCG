<<<<<<< HEAD
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

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
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(nonbonded)
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

# add target that generates version.c every time a make is run
# only do this if we generate the version
if(USE_VERSION_H)
    add_custom_target(gmx_version ALL
            COMMAND ${CMAKE_COMMAND} 
                -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                -D GIT_VERSION="${GIT_VERSION}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_C_CMAKEIN="${CMAKE_SOURCE_DIR}/src/gmxlib/version.c.cmakein"
                -D VERSION_C_OUT="${CMAKE_CURRENT_BINARY_DIR}/version.c"
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake 
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/gmxlib 
            DEPENDS ${CMAKE_SOURCE_DIR}/src/gmxlib/version.c.cmakein
            COMMENT "Generating version information")
endif()

# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.
<<<<<<< HEAD
file(GLOB GMXLIB_SOURCES *.c 
     selection/*.c trajana/*.c
     statistics/*.c nonbonded/*.c nonbonded/nb_kernel_c/*.c
     nonbonded/nb_kernel_adress_c/*.c)
=======
file(GLOB GMXLIB_SOURCES *.c selection/*.c trajana/*.c statistics/*.c)

if(NOT GMX_USE_PLUGINS)
  list(REMOVE_ITEM GMXLIB_SOURCES vmdio.c vmddlopen.c)
endif()
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

# This source file is generated
file(GLOB VERSION_SOURCE version.c)

if(VERSION_SOURCE)
	list(REMOVE_ITEM GMXLIB_SOURCES ${VERSION_SOURCE})
endif(VERSION_SOURCE)

# add version.c to the list of sources and tell cmake that it is generated
if(USE_VERSION_H)
LIST(APPEND GMXLIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/version.c) # auto-generated
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/version.c 
                PROPERTIES GENERATED true)
endif()

<<<<<<< HEAD
if(GMX_DOUBLE)
  set(SSETYPE sse2)
else()
  set(SSETYPE sse)
endif()

if(GMX_IA32_ASM)
  file(GLOB GMX_SSEKERNEL_C_SRC   nonbonded/nb_kernel_ia32_${SSETYPE}/*.c)
  if(GMX_ASM_USEASM-NASM)
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_ia32_${SSETYPE}/*intel_syntax*.s)    
  else()
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_ia32_${SSETYPE}/*${SSETYPE}.s nonbonded/nb_kernel_ia32_${SSETYPE}/*asm.s)
  endif()
endif(GMX_IA32_ASM)

if(GMX_X86_64_ASM)
  file(GLOB GMX_SSEKERNEL_C_SRC   nonbonded/nb_kernel_x86_64_${SSETYPE}/*.c)
  if(GMX_ASM_USEASM-NASM)
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_x86_64_${SSETYPE}/*intel_syntax*.s)
  else()
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_x86_64_${SSETYPE}/*${SSETYPE}.s nonbonded/nb_kernel_x86_64_${SSETYPE}/*asm.s)
  endif()
endif(GMX_X86_64_ASM)

if(GMX_FORTRAN)
  if (GMX_DOUBLE)
    file(GLOB FORTRAN_SOURCES nonbonded/nb_kernel_f77_double/*.[cf])
  else(GMX_DOUBLE)
    file(GLOB FORTRAN_SOURCES nonbonded/nb_kernel_f77_single/*.[cf])
  endif(GMX_DOUBLE)
endif(GMX_FORTRAN)

if(GMX_POWER6)
  file(GLOB FORTRAN_SOURCES nonbonded/nb_kernel_power6/*.[cF])
endif(GMX_POWER6)

if(GMX_BLUEGENE)
  file(GLOB GMX_BLUEGENE_C_SRC nonbonded/nb_kernel_bluegene/*.c)
endif(GMX_BLUEGENE)

if(GMX_PPC_ALTIVEC)
  file(GLOB GMX_PPC_ALTIVEC_SRC nonbonded/nb_kernel_ppc_altivec/*.c)
endif(GMX_PPC_ALTIVEC)

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
if(NOT GMX_EXTERNAL_BLAS)
  file(GLOB BLAS_SOURCES gmx_blas/*.c)
endif(NOT GMX_EXTERNAL_BLAS)

if(NOT GMX_EXTERNAL_LAPACK)
  file(GLOB LAPACK_SOURCES gmx_lapack/*.c)
endif(NOT GMX_EXTERNAL_LAPACK)

# This would be the standard way to include thread_mpi, but we want libgmx
# to link the functions directly
<<<<<<< HEAD
#if(GMX_THREADS)
#    add_subdirectory(thread_mpi)
#endif(GMX_THREADS)
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_MPI_LIB})
=======
#if(GMX_THREAD_MPI)
#    add_subdirectory(thread_mpi)
#endif(GMX_THREAD_MPI)
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${THREAD_MPI_LIB})

# apply gcc 4.4.x bug workaround
if(GMX_USE_GCC44_BUG_WORKAROUND)
   include(gmxGCC44O3BugWorkaround)
   gmx_apply_gcc44_bug_workaround("bondfree.c")
endif()
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

# Files called xxx_test.c are test drivers with a main() function for module xxx.c,
# so they should not be included in the library
file(GLOB_RECURSE NOT_GMXLIB_SOURCES *_test.c *\#*)
list(REMOVE_ITEM GMXLIB_SOURCES ${NOT_GMXLIB_SOURCES})  
# Selection has test_ instead of _test.
file(GLOB SELECTION_TEST selection/test*)
list(REMOVE_ITEM GMXLIB_SOURCES ${SELECTION_TEST})

<<<<<<< HEAD
# only fiddle with assembly kernels if we're not doing OpenMM build
if(NOT GMX_OPENMM) 
if(GMX_ASM_USEASM-NASM)
  enable_language(ASM-NASM)

  # if NASM is used, we need a special build command for windows...
  FOREACH(SRC ${GMX_SSEKERNEL_ASM_SRC})
    GET_FILENAME_COMPONENT(FILE_BASE ${SRC} NAME_WE)
    SET(OBJ ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}${CMAKE_C_OUTPUT_EXTENSION})

    ADD_CUSTOM_COMMAND(OUTPUT ${OBJ}
                       MAIN_DEPENDENCY ${SRC}
                       COMMAND ${CMAKE_ASM-NASM_COMPILER} -f ${CMAKE_ASM-NASM_OBJECT_FORMAT} -o ${OBJ} ${SRC})

    SET(ALL_ASM_OBJS ${ALL_ASM_OBJS} ${OBJ})
  ENDFOREACH(SRC ${GMX_SSEKERNEL_ASM_SRC})
  set(GMX_SSEKERNEL_ASM_SRC ${ALL_ASM_OBJS})

else(GMX_ASM_USEASM-NASM)

  enable_language(ASM-ATT)
  SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
  if(GMX_IA32_ASM)
    set_source_files_properties(${GMX_SSEKERNEL_ASM_SRC} PROPERTIES COMPILE_FLAGS "-c -m32")
  else()
    set_source_files_properties(${GMX_SSEKERNEL_ASM_SRC} PROPERTIES COMPILE_FLAGS "-c -m64")
  endif()

endif(GMX_ASM_USEASM-NASM)
endif(NOT GMX_OPENMM)

add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${GMX_SSEKERNEL_C_SRC} ${GMX_SSEKERNEL_ASM_SRC} ${FORTRAN_SOURCES} ${GMX_BLUEGENE_C_SRC} ${GMX_PPC_ALTIVEC_SRC} ${THREAD_MPI_SRC})
target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES}  ${THREAD_LIB})
if(USE_VERSION_H)
	add_dependencies(gmx gmx_version) 
endif()
set_target_properties(gmx PROPERTIES OUTPUT_NAME "gmx${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")
=======
# gpu utils + cuda tools module
if(GMX_GPU)
    add_subdirectory(cuda_tools)
    add_subdirectory(gpu_utils)   
    set(GMX_GPU_LIBRARIES gpu_utils cuda_tools)
endif()

# NONBONDED_SOURCES is imported from the nonbonded subdirectory.
add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${THREAD_MPI_SRC} ${NONBONDED_SOURCES})
target_link_libraries(gmx ${GMX_GPU_LIBRARIES} ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${THREAD_LIB} ${OpenMP_SHARED_LINKER_FLAGS})
if(USE_VERSION_H)
	add_dependencies(gmx gmx_version) 
endif()
set_target_properties(gmx PROPERTIES OUTPUT_NAME "gmx${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}"
    COMPILE_FLAGS "${OpenMP_C_FLAGS}" COMPILE_DEFINITIONS "TMPI_EXPORTS")
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmx.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmx${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)

