# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build

# Include any dependencies generated for this target.
include CMakeFiles/ClebschGordanTest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ClebschGordanTest.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ClebschGordanTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ClebschGordanTest.dir/flags.make

CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o: CMakeFiles/ClebschGordanTest.dir/flags.make
CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o: /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/Tests/ClebschGordanTest.cpp
CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o: CMakeFiles/ClebschGordanTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o -MF CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o.d -o CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o -c /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/Tests/ClebschGordanTest.cpp

CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/Tests/ClebschGordanTest.cpp > CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.i

CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/Tests/ClebschGordanTest.cpp -o CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.s

# Object files for target ClebschGordanTest
ClebschGordanTest_OBJECTS = \
"CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o"

# External object files for target ClebschGordanTest
ClebschGordanTest_EXTERNAL_OBJECTS =

ClebschGordanTest: CMakeFiles/ClebschGordanTest.dir/Tests/ClebschGordanTest.cpp.o
ClebschGordanTest: CMakeFiles/ClebschGordanTest.dir/build.make
ClebschGordanTest: /usr/lib/x86_64-linux-gnu/libgtest.a
ClebschGordanTest: libThreeBodyDecaysLib.a
ClebschGordanTest: CMakeFiles/ClebschGordanTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ClebschGordanTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ClebschGordanTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ClebschGordanTest.dir/build: ClebschGordanTest
.PHONY : CMakeFiles/ClebschGordanTest.dir/build

CMakeFiles/ClebschGordanTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ClebschGordanTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ClebschGordanTest.dir/clean

CMakeFiles/ClebschGordanTest.dir/depend:
	cd /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build /home/hvdsmagt/ThreeBodyDecays.jl/ThreeBodyDecays/build/CMakeFiles/ClebschGordanTest.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/ClebschGordanTest.dir/depend

