# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.1

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
CMAKE_COMMAND = /panfs-hydra/user/home/mgrecu/cmake-3.1.0-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /panfs-hydra/user/home/mgrecu/cmake-3.1.0-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /PANFS/user/home/mgrecu/data/flann-1.8.4-src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/flann_example_c.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/flann_example_c.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/flann_example_c.dir/flags.make

examples/CMakeFiles/flann_example_c.dir/flann_example.c.o: examples/CMakeFiles/flann_example_c.dir/flags.make
examples/CMakeFiles/flann_example_c.dir/flann_example.c.o: ../examples/flann_example.c
	$(CMAKE_COMMAND) -E cmake_progress_report /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object examples/CMakeFiles/flann_example_c.dir/flann_example.c.o"
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/flann_example_c.dir/flann_example.c.o   -c /PANFS/user/home/mgrecu/data/flann-1.8.4-src/examples/flann_example.c

examples/CMakeFiles/flann_example_c.dir/flann_example.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/flann_example_c.dir/flann_example.c.i"
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /PANFS/user/home/mgrecu/data/flann-1.8.4-src/examples/flann_example.c > CMakeFiles/flann_example_c.dir/flann_example.c.i

examples/CMakeFiles/flann_example_c.dir/flann_example.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/flann_example_c.dir/flann_example.c.s"
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /PANFS/user/home/mgrecu/data/flann-1.8.4-src/examples/flann_example.c -o CMakeFiles/flann_example_c.dir/flann_example.c.s

examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.requires:
.PHONY : examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.requires

examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.provides: examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.requires
	$(MAKE) -f examples/CMakeFiles/flann_example_c.dir/build.make examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.provides.build
.PHONY : examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.provides

examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.provides.build: examples/CMakeFiles/flann_example_c.dir/flann_example.c.o

# Object files for target flann_example_c
flann_example_c_OBJECTS = \
"CMakeFiles/flann_example_c.dir/flann_example.c.o"

# External object files for target flann_example_c
flann_example_c_EXTERNAL_OBJECTS =

bin/flann_example_c: examples/CMakeFiles/flann_example_c.dir/flann_example.c.o
bin/flann_example_c: examples/CMakeFiles/flann_example_c.dir/build.make
bin/flann_example_c: lib/libflann.so.1.8.4
bin/flann_example_c: lib/libflann_s.a
bin/flann_example_c: examples/CMakeFiles/flann_example_c.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/flann_example_c"
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/flann_example_c.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/flann_example_c.dir/build: bin/flann_example_c
.PHONY : examples/CMakeFiles/flann_example_c.dir/build

examples/CMakeFiles/flann_example_c.dir/requires: examples/CMakeFiles/flann_example_c.dir/flann_example.c.o.requires
.PHONY : examples/CMakeFiles/flann_example_c.dir/requires

examples/CMakeFiles/flann_example_c.dir/clean:
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/flann_example_c.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/flann_example_c.dir/clean

examples/CMakeFiles/flann_example_c.dir/depend:
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /PANFS/user/home/mgrecu/data/flann-1.8.4-src /PANFS/user/home/mgrecu/data/flann-1.8.4-src/examples /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/examples/CMakeFiles/flann_example_c.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/flann_example_c.dir/depend

