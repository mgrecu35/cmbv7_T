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

# Utility rule file for auxclean.

# Include the progress variables for this target.
include doc/CMakeFiles/auxclean.dir/progress.make

doc/CMakeFiles/auxclean:
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc && /panfs-hydra/user/home/mgrecu/cmake-3.1.0-Linux-x86_64/bin/cmake -E remove /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc/manual.aux /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc/manual.idx /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc/manual.ind

auxclean: doc/CMakeFiles/auxclean
auxclean: doc/CMakeFiles/auxclean.dir/build.make
.PHONY : auxclean

# Rule to build all files generated by this target.
doc/CMakeFiles/auxclean.dir/build: auxclean
.PHONY : doc/CMakeFiles/auxclean.dir/build

doc/CMakeFiles/auxclean.dir/clean:
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc && $(CMAKE_COMMAND) -P CMakeFiles/auxclean.dir/cmake_clean.cmake
.PHONY : doc/CMakeFiles/auxclean.dir/clean

doc/CMakeFiles/auxclean.dir/depend:
	cd /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /PANFS/user/home/mgrecu/data/flann-1.8.4-src /PANFS/user/home/mgrecu/data/flann-1.8.4-src/doc /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc /PANFS/user/home/mgrecu/data/flann-1.8.4-src/build/doc/CMakeFiles/auxclean.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/CMakeFiles/auxclean.dir/depend

