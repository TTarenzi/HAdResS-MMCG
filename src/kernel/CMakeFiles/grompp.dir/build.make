# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/pckr144/potestio/hadressmacs/adressmacs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/pckr144/potestio/hadressmacs/adressmacs

# Include any dependencies generated for this target.
include src/kernel/CMakeFiles/grompp.dir/depend.make

# Include the progress variables for this target.
include src/kernel/CMakeFiles/grompp.dir/progress.make

# Include the compile flags for this target's objects.
include src/kernel/CMakeFiles/grompp.dir/flags.make

src/kernel/CMakeFiles/grompp.dir/grompp.c.o: src/kernel/CMakeFiles/grompp.dir/flags.make
src/kernel/CMakeFiles/grompp.dir/grompp.c.o: src/kernel/grompp.c
	$(CMAKE_COMMAND) -E cmake_progress_report /data/pckr144/potestio/hadressmacs/adressmacs/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/kernel/CMakeFiles/grompp.dir/grompp.c.o"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/grompp.dir/grompp.c.o   -c /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel/grompp.c

src/kernel/CMakeFiles/grompp.dir/grompp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/grompp.dir/grompp.c.i"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel/grompp.c > CMakeFiles/grompp.dir/grompp.c.i

src/kernel/CMakeFiles/grompp.dir/grompp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/grompp.dir/grompp.c.s"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel/grompp.c -o CMakeFiles/grompp.dir/grompp.c.s

src/kernel/CMakeFiles/grompp.dir/grompp.c.o.requires:
.PHONY : src/kernel/CMakeFiles/grompp.dir/grompp.c.o.requires

src/kernel/CMakeFiles/grompp.dir/grompp.c.o.provides: src/kernel/CMakeFiles/grompp.dir/grompp.c.o.requires
	$(MAKE) -f src/kernel/CMakeFiles/grompp.dir/build.make src/kernel/CMakeFiles/grompp.dir/grompp.c.o.provides.build
.PHONY : src/kernel/CMakeFiles/grompp.dir/grompp.c.o.provides

src/kernel/CMakeFiles/grompp.dir/grompp.c.o.provides.build: src/kernel/CMakeFiles/grompp.dir/grompp.c.o

# Object files for target grompp
grompp_OBJECTS = \
"CMakeFiles/grompp.dir/grompp.c.o"

# External object files for target grompp
grompp_EXTERNAL_OBJECTS =

src/kernel/grompp: src/kernel/CMakeFiles/grompp.dir/grompp.c.o
src/kernel/grompp: src/kernel/CMakeFiles/grompp.dir/build.make
src/kernel/grompp: src/kernel/libgmxpreprocess.so.6
src/kernel/grompp: src/mdlib/libmd.so.6
src/kernel/grompp: src/gmxlib/libgmx.so.6
src/kernel/grompp: /usr/lib64/libfftw3f.so
src/kernel/grompp: /usr/lib64/libxml2.so
src/kernel/grompp: src/kernel/CMakeFiles/grompp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable grompp"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/grompp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/kernel/CMakeFiles/grompp.dir/build: src/kernel/grompp
.PHONY : src/kernel/CMakeFiles/grompp.dir/build

# Object files for target grompp
grompp_OBJECTS = \
"CMakeFiles/grompp.dir/grompp.c.o"

# External object files for target grompp
grompp_EXTERNAL_OBJECTS =

src/kernel/CMakeFiles/CMakeRelink.dir/grompp: src/kernel/CMakeFiles/grompp.dir/grompp.c.o
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: src/kernel/CMakeFiles/grompp.dir/build.make
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: src/kernel/libgmxpreprocess.so.6
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: src/mdlib/libmd.so.6
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: src/gmxlib/libgmx.so.6
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: /usr/lib64/libfftw3f.so
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: /usr/lib64/libxml2.so
src/kernel/CMakeFiles/CMakeRelink.dir/grompp: src/kernel/CMakeFiles/grompp.dir/relink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable CMakeFiles/CMakeRelink.dir/grompp"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/grompp.dir/relink.txt --verbose=$(VERBOSE)

# Rule to relink during preinstall.
src/kernel/CMakeFiles/grompp.dir/preinstall: src/kernel/CMakeFiles/CMakeRelink.dir/grompp
.PHONY : src/kernel/CMakeFiles/grompp.dir/preinstall

src/kernel/CMakeFiles/grompp.dir/requires: src/kernel/CMakeFiles/grompp.dir/grompp.c.o.requires
.PHONY : src/kernel/CMakeFiles/grompp.dir/requires

src/kernel/CMakeFiles/grompp.dir/clean:
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel && $(CMAKE_COMMAND) -P CMakeFiles/grompp.dir/cmake_clean.cmake
.PHONY : src/kernel/CMakeFiles/grompp.dir/clean

src/kernel/CMakeFiles/grompp.dir/depend:
	cd /data/pckr144/potestio/hadressmacs/adressmacs && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/pckr144/potestio/hadressmacs/adressmacs /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel /data/pckr144/potestio/hadressmacs/adressmacs /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel /data/pckr144/potestio/hadressmacs/adressmacs/src/kernel/CMakeFiles/grompp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/kernel/CMakeFiles/grompp.dir/depend

