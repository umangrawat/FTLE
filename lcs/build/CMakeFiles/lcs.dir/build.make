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
CMAKE_SOURCE_DIR = /export/home/urawat/Desktop/lcs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /export/home/urawat/Desktop/lcs/build

# Include any dependencies generated for this target.
include CMakeFiles/lcs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lcs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lcs.dir/flags.make

CMakeFiles/lcs.dir/lcs.cpp.o: CMakeFiles/lcs.dir/flags.make
CMakeFiles/lcs.dir/lcs.cpp.o: ../lcs.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/export/home/urawat/Desktop/lcs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lcs.dir/lcs.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lcs.dir/lcs.cpp.o -c /export/home/urawat/Desktop/lcs/lcs.cpp

CMakeFiles/lcs.dir/lcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lcs.dir/lcs.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /export/home/urawat/Desktop/lcs/lcs.cpp > CMakeFiles/lcs.dir/lcs.cpp.i

CMakeFiles/lcs.dir/lcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lcs.dir/lcs.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /export/home/urawat/Desktop/lcs/lcs.cpp -o CMakeFiles/lcs.dir/lcs.cpp.s

CMakeFiles/lcs.dir/lcs.cpp.o.requires:

.PHONY : CMakeFiles/lcs.dir/lcs.cpp.o.requires

CMakeFiles/lcs.dir/lcs.cpp.o.provides: CMakeFiles/lcs.dir/lcs.cpp.o.requires
	$(MAKE) -f CMakeFiles/lcs.dir/build.make CMakeFiles/lcs.dir/lcs.cpp.o.provides.build
.PHONY : CMakeFiles/lcs.dir/lcs.cpp.o.provides

CMakeFiles/lcs.dir/lcs.cpp.o.provides.build: CMakeFiles/lcs.dir/lcs.cpp.o


# Object files for target lcs
lcs_OBJECTS = \
"CMakeFiles/lcs.dir/lcs.cpp.o"

# External object files for target lcs
lcs_EXTERNAL_OBJECTS =

lcs: CMakeFiles/lcs.dir/lcs.cpp.o
lcs: CMakeFiles/lcs.dir/build.make
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOTecplotTable-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersSelection-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOExodus-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkLocalExample-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOVideo-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOImport-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkTestingGenericBridge-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersPoints-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOMINC-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOExportOpenGL2-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOExport-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingGL2PSOpenGL2-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkgl2ps-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersFlowPaths-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersProgrammable-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOInfovis-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtklibxml2-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkDomainsChemistryOpenGL2-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkDomainsChemistry-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingContextOpenGL2-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersSMP-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingImage-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingStencil-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersTexture-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersParallelImaging-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOAMR-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOEnSight-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersVerdict-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOPLY-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingMorphological-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkViewsContext2D-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkTestingIOSQL-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOSQL-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkInteractionImage-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOMovie-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkoggtheora-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingStatistics-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingLOD-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersHyperTree-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOLSDyna-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersGeneric-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOParallel-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkTestingRendering-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOParallelXML-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkGeovisCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersTopology-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkViewsInfovis-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingVolumeOpenGL2-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtklibharu-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersAMR-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkverdict-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtksqlite-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkexoIIc-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIONetCDF-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtknetcdfcpp-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkNetCDF-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkhdf5_hl-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkhdf5-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersParallel-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkjsoncpp-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOGeometry-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkParallelCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOLegacy-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkproj4-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkChartsCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingContext2D-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersImaging-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkViewsCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkInteractionWidgets-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingGeneral-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersHybrid-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingSources-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkInteractionStyle-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingAnnotation-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingColor-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingLabel-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingFreeType-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkfreetype-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkInfovisLayout-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingHybrid-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOImage-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkDICOMParser-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkmetaio-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkpng-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtktiff-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkjpeg-8.1.so.1
lcs: /usr/lib/x86_64-linux-gnu/libm.so
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersModeling-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkInfovisCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersExtraction-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersStatistics-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingFourier-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkalglib-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingOpenGL2-8.1.so.1
lcs: /usr/lib/x86_64-linux-gnu/libSM.so
lcs: /usr/lib/x86_64-linux-gnu/libICE.so
lcs: /usr/lib/x86_64-linux-gnu/libX11.so
lcs: /usr/lib/x86_64-linux-gnu/libXext.so
lcs: /usr/lib/x86_64-linux-gnu/libXt.so
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkglew-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingVolume-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkRenderingCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonColor-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersGeometry-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersSources-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersGeneral-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonComputationalGeometry-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkFiltersCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOXML-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOXMLParser-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkIOCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkzlib-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtklz4-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkexpat-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkImagingMath-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonExecutionModel-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonDataModel-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonTransforms-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonMisc-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonMath-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonSystem-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtkCommonCore-8.1.so.1
lcs: /export/home/urawat/projects/VTK-build/lib/libvtksys-8.1.so.1
lcs: CMakeFiles/lcs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/export/home/urawat/Desktop/lcs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lcs"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lcs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lcs.dir/build: lcs

.PHONY : CMakeFiles/lcs.dir/build

CMakeFiles/lcs.dir/requires: CMakeFiles/lcs.dir/lcs.cpp.o.requires

.PHONY : CMakeFiles/lcs.dir/requires

CMakeFiles/lcs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lcs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lcs.dir/clean

CMakeFiles/lcs.dir/depend:
	cd /export/home/urawat/Desktop/lcs/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /export/home/urawat/Desktop/lcs /export/home/urawat/Desktop/lcs /export/home/urawat/Desktop/lcs/build /export/home/urawat/Desktop/lcs/build /export/home/urawat/Desktop/lcs/build/CMakeFiles/lcs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lcs.dir/depend

