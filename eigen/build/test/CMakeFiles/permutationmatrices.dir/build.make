# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /export/home/urawat/Desktop/eigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /export/home/urawat/Desktop/eigen/build

# Utility rule file for permutationmatrices.

# Include the progress variables for this target.
include test/CMakeFiles/permutationmatrices.dir/progress.make

permutationmatrices: test/CMakeFiles/permutationmatrices.dir/build.make

.PHONY : permutationmatrices

# Rule to build all files generated by this target.
test/CMakeFiles/permutationmatrices.dir/build: permutationmatrices

.PHONY : test/CMakeFiles/permutationmatrices.dir/build

test/CMakeFiles/permutationmatrices.dir/clean:
	cd /export/home/urawat/Desktop/eigen/build/test && $(CMAKE_COMMAND) -P CMakeFiles/permutationmatrices.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/permutationmatrices.dir/clean

test/CMakeFiles/permutationmatrices.dir/depend:
	cd /export/home/urawat/Desktop/eigen/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /export/home/urawat/Desktop/eigen /export/home/urawat/Desktop/eigen/test /export/home/urawat/Desktop/eigen/build /export/home/urawat/Desktop/eigen/build/test /export/home/urawat/Desktop/eigen/build/test/CMakeFiles/permutationmatrices.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/permutationmatrices.dir/depend

