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
include src/tools/CMakeFiles/g_bond.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/g_bond.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/g_bond.dir/flags.make

src/tools/CMakeFiles/g_bond.dir/g_bond.c.o: src/tools/CMakeFiles/g_bond.dir/flags.make
src/tools/CMakeFiles/g_bond.dir/g_bond.c.o: src/tools/g_bond.c
	$(CMAKE_COMMAND) -E cmake_progress_report /data/pckr144/potestio/hadressmacs/adressmacs/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/tools/CMakeFiles/g_bond.dir/g_bond.c.o"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/g_bond.dir/g_bond.c.o   -c /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/g_bond.c

src/tools/CMakeFiles/g_bond.dir/g_bond.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/g_bond.dir/g_bond.c.i"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/g_bond.c > CMakeFiles/g_bond.dir/g_bond.c.i

src/tools/CMakeFiles/g_bond.dir/g_bond.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/g_bond.dir/g_bond.c.s"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/g_bond.c -o CMakeFiles/g_bond.dir/g_bond.c.s

src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.requires:
.PHONY : src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.requires

src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.provides: src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.requires
	$(MAKE) -f src/tools/CMakeFiles/g_bond.dir/build.make src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.provides.build
.PHONY : src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.provides

src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.provides.build: src/tools/CMakeFiles/g_bond.dir/g_bond.c.o

# Object files for target g_bond
g_bond_OBJECTS = \
"CMakeFiles/g_bond.dir/g_bond.c.o"

# External object files for target g_bond
g_bond_EXTERNAL_OBJECTS =

src/tools/g_bond: src/tools/CMakeFiles/g_bond.dir/g_bond.c.o
src/tools/g_bond: src/tools/CMakeFiles/g_bond.dir/build.make
src/tools/g_bond: src/tools/libgmxana.so.6
src/tools/g_bond: src/mdlib/libmd.so.6
src/tools/g_bond: /usr/lib64/libfftw3f.so
src/tools/g_bond: /usr/lib64/libxml2.so
src/tools/g_bond: src/gmxlib/libgmx.so.6
src/tools/g_bond: src/tools/CMakeFiles/g_bond.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable g_bond"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/g_bond.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/g_bond.dir/build: src/tools/g_bond
.PHONY : src/tools/CMakeFiles/g_bond.dir/build

# Object files for target g_bond
g_bond_OBJECTS = \
"CMakeFiles/g_bond.dir/g_bond.c.o"

# External object files for target g_bond
g_bond_EXTERNAL_OBJECTS =

src/tools/CMakeFiles/CMakeRelink.dir/g_bond: src/tools/CMakeFiles/g_bond.dir/g_bond.c.o
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: src/tools/CMakeFiles/g_bond.dir/build.make
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: src/tools/libgmxana.so.6
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: src/mdlib/libmd.so.6
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: /usr/lib64/libfftw3f.so
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: /usr/lib64/libxml2.so
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: src/gmxlib/libgmx.so.6
src/tools/CMakeFiles/CMakeRelink.dir/g_bond: src/tools/CMakeFiles/g_bond.dir/relink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable CMakeFiles/CMakeRelink.dir/g_bond"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/g_bond.dir/relink.txt --verbose=$(VERBOSE)

# Rule to relink during preinstall.
src/tools/CMakeFiles/g_bond.dir/preinstall: src/tools/CMakeFiles/CMakeRelink.dir/g_bond
.PHONY : src/tools/CMakeFiles/g_bond.dir/preinstall

src/tools/CMakeFiles/g_bond.dir/requires: src/tools/CMakeFiles/g_bond.dir/g_bond.c.o.requires
.PHONY : src/tools/CMakeFiles/g_bond.dir/requires

src/tools/CMakeFiles/g_bond.dir/clean:
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/g_bond.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/g_bond.dir/clean

src/tools/CMakeFiles/g_bond.dir/depend:
	cd /data/pckr144/potestio/hadressmacs/adressmacs && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/pckr144/potestio/hadressmacs/adressmacs /data/pckr144/potestio/hadressmacs/adressmacs/src/tools /data/pckr144/potestio/hadressmacs/adressmacs /data/pckr144/potestio/hadressmacs/adressmacs/src/tools /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/CMakeFiles/g_bond.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/g_bond.dir/depend

