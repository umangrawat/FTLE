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

# Include any dependencies generated for this target.
include test/CMakeFiles/sparselu_3.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/sparselu_3.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/sparselu_3.dir/flags.make

test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o: test/CMakeFiles/sparselu_3.dir/flags.make
test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o: ../test/sparselu.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/export/home/urawat/Desktop/eigen/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o"
	cd /export/home/urawat/Desktop/eigen/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sparselu_3.dir/sparselu.cpp.o -c /export/home/urawat/Desktop/eigen/test/sparselu.cpp

test/CMakeFiles/sparselu_3.dir/sparselu.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sparselu_3.dir/sparselu.cpp.i"
	cd /export/home/urawat/Desktop/eigen/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /export/home/urawat/Desktop/eigen/test/sparselu.cpp > CMakeFiles/sparselu_3.dir/sparselu.cpp.i

test/CMakeFiles/sparselu_3.dir/sparselu.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sparselu_3.dir/sparselu.cpp.s"
	cd /export/home/urawat/Desktop/eigen/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /export/home/urawat/Desktop/eigen/test/sparselu.cpp -o CMakeFiles/sparselu_3.dir/sparselu.cpp.s

test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.requires:

.PHONY : test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.requires

test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.provides: test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/sparselu_3.dir/build.make test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.provides.build
.PHONY : test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.provides

test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.provides.build: test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o


# Object files for target sparselu_3
sparselu_3_OBJECTS = \
"CMakeFiles/sparselu_3.dir/sparselu.cpp.o"

# External object files for target sparselu_3
sparselu_3_EXTERNAL_OBJECTS =

test/sparselu_3: test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o
test/sparselu_3: test/CMakeFiles/sparselu_3.dir/build.make
test/sparselu_3: test/CMakeFiles/sparselu_3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/export/home/urawat/Desktop/eigen/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sparselu_3"
	cd /export/home/urawat/Desktop/eigen/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sparselu_3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/sparselu_3.dir/build: test/sparselu_3

.PHONY : test/CMakeFiles/sparselu_3.dir/build

test/CMakeFiles/sparselu_3.dir/requires: test/CMakeFiles/sparselu_3.dir/sparselu.cpp.o.requires

.PHONY : test/CMakeFiles/sparselu_3.dir/requires

test/CMakeFiles/sparselu_3.dir/clean:
	cd /export/home/urawat/Desktop/eigen/build/test && $(CMAKE_COMMAND) -P CMakeFiles/sparselu_3.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/sparselu_3.dir/clean

test/CMakeFiles/sparselu_3.dir/depend:
	cd /export/home/urawat/Desktop/eigen/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /export/home/urawat/Desktop/eigen /export/home/urawat/Desktop/eigen/test /export/home/urawat/Desktop/eigen/build /export/home/urawat/Desktop/eigen/build/test /export/home/urawat/Desktop/eigen/build/test/CMakeFiles/sparselu_3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/sparselu_3.dir/depend

