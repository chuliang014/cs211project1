# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/assignmentCode.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/assignmentCode.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assignmentCode.dir/flags.make

CMakeFiles/assignmentCode.dir/main.c.o: CMakeFiles/assignmentCode.dir/flags.make
CMakeFiles/assignmentCode.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/assignmentCode.dir/main.c.o"
	/Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/assignmentCode.dir/main.c.o   -c "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/main.c"

CMakeFiles/assignmentCode.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/assignmentCode.dir/main.c.i"
	/Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/main.c" > CMakeFiles/assignmentCode.dir/main.c.i

CMakeFiles/assignmentCode.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/assignmentCode.dir/main.c.s"
	/Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/main.c" -o CMakeFiles/assignmentCode.dir/main.c.s

# Object files for target assignmentCode
assignmentCode_OBJECTS = \
"CMakeFiles/assignmentCode.dir/main.c.o"

# External object files for target assignmentCode
assignmentCode_EXTERNAL_OBJECTS =

assignmentCode: CMakeFiles/assignmentCode.dir/main.c.o
assignmentCode: CMakeFiles/assignmentCode.dir/build.make
assignmentCode: CMakeFiles/assignmentCode.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable assignmentCode"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assignmentCode.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assignmentCode.dir/build: assignmentCode

.PHONY : CMakeFiles/assignmentCode.dir/build

CMakeFiles/assignmentCode.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assignmentCode.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assignmentCode.dir/clean

CMakeFiles/assignmentCode.dir/depend:
	cd "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode" "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode" "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug" "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug" "/Users/vincent/Learning Data/cs211HighPerformanceComputing/assignmentCode/cmake-build-debug/CMakeFiles/assignmentCode.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/assignmentCode.dir/depend
