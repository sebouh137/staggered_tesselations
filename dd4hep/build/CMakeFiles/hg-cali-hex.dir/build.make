# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/cmake-3.26.3-2qt5tkelo3cojdo56hfe5cai2llcqxkj/bin/cmake

# The command to remove a file.
RM = /opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/cmake-3.26.3-2qt5tkelo3cojdo56hfe5cai2llcqxkj/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sebouh/staggered_tesselations/dd4hep

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sebouh/staggered_tesselations/dd4hep/build

# Include any dependencies generated for this target.
include CMakeFiles/hg-cali-hex.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/hg-cali-hex.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/hg-cali-hex.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/hg-cali-hex.dir/flags.make

CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o: CMakeFiles/hg-cali-hex.dir/flags.make
CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o: /home/sebouh/staggered_tesselations/dd4hep/src/HexGrid.cpp
CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o: CMakeFiles/hg-cali-hex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sebouh/staggered_tesselations/dd4hep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o -MF CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o.d -o CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o -c /home/sebouh/staggered_tesselations/dd4hep/src/HexGrid.cpp

CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sebouh/staggered_tesselations/dd4hep/src/HexGrid.cpp > CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.i

CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sebouh/staggered_tesselations/dd4hep/src/HexGrid.cpp -o CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.s

CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o: CMakeFiles/hg-cali-hex.dir/flags.make
CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o: /home/sebouh/staggered_tesselations/dd4hep/src/InsertCalorimeter_geo.cpp
CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o: CMakeFiles/hg-cali-hex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sebouh/staggered_tesselations/dd4hep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o -MF CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o.d -o CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o -c /home/sebouh/staggered_tesselations/dd4hep/src/InsertCalorimeter_geo.cpp

CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sebouh/staggered_tesselations/dd4hep/src/InsertCalorimeter_geo.cpp > CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.i

CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sebouh/staggered_tesselations/dd4hep/src/InsertCalorimeter_geo.cpp -o CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.s

# Object files for target hg-cali-hex
hg__cali__hex_OBJECTS = \
"CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o" \
"CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o"

# External object files for target hg-cali-hex
hg__cali__hex_EXTERNAL_OBJECTS =

libhg-cali-hex.so: CMakeFiles/hg-cali-hex.dir/src/HexGrid.cpp.o
libhg-cali-hex.so: CMakeFiles/hg-cali-hex.dir/src/InsertCalorimeter_geo.cpp.o
libhg-cali-hex.so: CMakeFiles/hg-cali-hex.dir/build.make
libhg-cali-hex.so: /usr/local/lib/libDDRec.so.1.25
libhg-cali-hex.so: /usr/local/lib/libDDCore.so.1.25
libhg-cali-hex.so: /usr/local/lib/libDDParsers.so.1.25
libhg-cali-hex.so: /usr/local/lib/root/libRint.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libTree.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libPhysics.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libGeom.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libHist.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libMatrix.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libGenVector.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libMathCore.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libImt.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libMultiProc.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libNet.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libRIO.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libThread.so.6.26.10
libhg-cali-hex.so: /usr/local/lib/root/libCore.so.6.26.10
libhg-cali-hex.so: /opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/xerces-c-3.2.4-g233npbtqnly6pxggbwrsuplsamqbyqb/lib/libxerces-c.so
libhg-cali-hex.so: CMakeFiles/hg-cali-hex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sebouh/staggered_tesselations/dd4hep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library libhg-cali-hex.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hg-cali-hex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/hg-cali-hex.dir/build: libhg-cali-hex.so
.PHONY : CMakeFiles/hg-cali-hex.dir/build

CMakeFiles/hg-cali-hex.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/hg-cali-hex.dir/cmake_clean.cmake
.PHONY : CMakeFiles/hg-cali-hex.dir/clean

CMakeFiles/hg-cali-hex.dir/depend:
	cd /home/sebouh/staggered_tesselations/dd4hep/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sebouh/staggered_tesselations/dd4hep /home/sebouh/staggered_tesselations/dd4hep /home/sebouh/staggered_tesselations/dd4hep/build /home/sebouh/staggered_tesselations/dd4hep/build /home/sebouh/staggered_tesselations/dd4hep/build/CMakeFiles/hg-cali-hex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/hg-cali-hex.dir/depend

