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
include src/tools/CMakeFiles/g_tcaf.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/g_tcaf.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/g_tcaf.dir/flags.make

src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o: src/tools/CMakeFiles/g_tcaf.dir/flags.make
src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o: src/tools/g_tcaf.c
	$(CMAKE_COMMAND) -E cmake_progress_report /data/pckr144/potestio/hadressmacs/adressmacs/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/g_tcaf.dir/g_tcaf.c.o   -c /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/g_tcaf.c

src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/g_tcaf.dir/g_tcaf.c.i"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/g_tcaf.c > CMakeFiles/g_tcaf.dir/g_tcaf.c.i

src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/g_tcaf.dir/g_tcaf.c.s"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/g_tcaf.c -o CMakeFiles/g_tcaf.dir/g_tcaf.c.s

src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.requires:
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.requires

src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.provides: src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.requires
	$(MAKE) -f src/tools/CMakeFiles/g_tcaf.dir/build.make src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.provides.build
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.provides

src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.provides.build: src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o

# Object files for target g_tcaf
g_tcaf_OBJECTS = \
"CMakeFiles/g_tcaf.dir/g_tcaf.c.o"

# External object files for target g_tcaf
g_tcaf_EXTERNAL_OBJECTS =

src/tools/g_tcaf: src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o
src/tools/g_tcaf: src/tools/CMakeFiles/g_tcaf.dir/build.make
src/tools/g_tcaf: src/tools/libgmxana.so.6
src/tools/g_tcaf: src/mdlib/libmd.so.6
src/tools/g_tcaf: /usr/lib64/libfftw3f.so
src/tools/g_tcaf: /usr/lib64/libxml2.so
src/tools/g_tcaf: src/gmxlib/libgmx.so.6
src/tools/g_tcaf: src/tools/CMakeFiles/g_tcaf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable g_tcaf"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/g_tcaf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/g_tcaf.dir/build: src/tools/g_tcaf
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/build

# Object files for target g_tcaf
g_tcaf_OBJECTS = \
"CMakeFiles/g_tcaf.dir/g_tcaf.c.o"

# External object files for target g_tcaf
g_tcaf_EXTERNAL_OBJECTS =

src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: src/tools/CMakeFiles/g_tcaf.dir/build.make
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: src/tools/libgmxana.so.6
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: src/mdlib/libmd.so.6
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: /usr/lib64/libfftw3f.so
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: /usr/lib64/libxml2.so
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: src/gmxlib/libgmx.so.6
src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf: src/tools/CMakeFiles/g_tcaf.dir/relink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable CMakeFiles/CMakeRelink.dir/g_tcaf"
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/g_tcaf.dir/relink.txt --verbose=$(VERBOSE)

# Rule to relink during preinstall.
src/tools/CMakeFiles/g_tcaf.dir/preinstall: src/tools/CMakeFiles/CMakeRelink.dir/g_tcaf
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/preinstall

src/tools/CMakeFiles/g_tcaf.dir/requires: src/tools/CMakeFiles/g_tcaf.dir/g_tcaf.c.o.requires
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/requires

src/tools/CMakeFiles/g_tcaf.dir/clean:
	cd /data/pckr144/potestio/hadressmacs/adressmacs/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/g_tcaf.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/clean

src/tools/CMakeFiles/g_tcaf.dir/depend:
	cd /data/pckr144/potestio/hadressmacs/adressmacs && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/pckr144/potestio/hadressmacs/adressmacs /data/pckr144/potestio/hadressmacs/adressmacs/src/tools /data/pckr144/potestio/hadressmacs/adressmacs /data/pckr144/potestio/hadressmacs/adressmacs/src/tools /data/pckr144/potestio/hadressmacs/adressmacs/src/tools/CMakeFiles/g_tcaf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/g_tcaf.dir/depend

