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
CMAKE_COMMAND = /home/yanglaoyuan/.local/lib/python3.8/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/yanglaoyuan/.local/lib/python3.8/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/d/Github/GraphPi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/Github/GraphPi

# Include any dependencies generated for this target.
include test/CMakeFiles/asgs_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/asgs_test.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/asgs_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/asgs_test.dir/flags.make

test/CMakeFiles/asgs_test.dir/gtest_main.cpp.o: test/CMakeFiles/asgs_test.dir/flags.make
test/CMakeFiles/asgs_test.dir/gtest_main.cpp.o: test/gtest_main.cpp
test/CMakeFiles/asgs_test.dir/gtest_main.cpp.o: test/CMakeFiles/asgs_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Github/GraphPi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/asgs_test.dir/gtest_main.cpp.o"
	cd /mnt/d/Github/GraphPi/test && mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/asgs_test.dir/gtest_main.cpp.o -MF CMakeFiles/asgs_test.dir/gtest_main.cpp.o.d -o CMakeFiles/asgs_test.dir/gtest_main.cpp.o -c /mnt/d/Github/GraphPi/test/gtest_main.cpp

test/CMakeFiles/asgs_test.dir/gtest_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/asgs_test.dir/gtest_main.cpp.i"
	cd /mnt/d/Github/GraphPi/test && mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Github/GraphPi/test/gtest_main.cpp > CMakeFiles/asgs_test.dir/gtest_main.cpp.i

test/CMakeFiles/asgs_test.dir/gtest_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/asgs_test.dir/gtest_main.cpp.s"
	cd /mnt/d/Github/GraphPi/test && mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Github/GraphPi/test/gtest_main.cpp -o CMakeFiles/asgs_test.dir/gtest_main.cpp.s

test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.o: test/CMakeFiles/asgs_test.dir/flags.make
test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.o: test/asgs_testing.cpp
test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.o: test/CMakeFiles/asgs_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Github/GraphPi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.o"
	cd /mnt/d/Github/GraphPi/test && mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.o -MF CMakeFiles/asgs_test.dir/asgs_testing.cpp.o.d -o CMakeFiles/asgs_test.dir/asgs_testing.cpp.o -c /mnt/d/Github/GraphPi/test/asgs_testing.cpp

test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/asgs_test.dir/asgs_testing.cpp.i"
	cd /mnt/d/Github/GraphPi/test && mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Github/GraphPi/test/asgs_testing.cpp > CMakeFiles/asgs_test.dir/asgs_testing.cpp.i

test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/asgs_test.dir/asgs_testing.cpp.s"
	cd /mnt/d/Github/GraphPi/test && mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Github/GraphPi/test/asgs_testing.cpp -o CMakeFiles/asgs_test.dir/asgs_testing.cpp.s

# Object files for target asgs_test
asgs_test_OBJECTS = \
"CMakeFiles/asgs_test.dir/gtest_main.cpp.o" \
"CMakeFiles/asgs_test.dir/asgs_testing.cpp.o"

# External object files for target asgs_test
asgs_test_EXTERNAL_OBJECTS =

bin/asgs_test: test/CMakeFiles/asgs_test.dir/gtest_main.cpp.o
bin/asgs_test: test/CMakeFiles/asgs_test.dir/asgs_testing.cpp.o
bin/asgs_test: test/CMakeFiles/asgs_test.dir/build.make
bin/asgs_test: /usr/local/lib/libgtest.a
bin/asgs_test: libs/libgraph_mining.so
bin/asgs_test: test/CMakeFiles/asgs_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/d/Github/GraphPi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../bin/asgs_test"
	cd /mnt/d/Github/GraphPi/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/asgs_test.dir/link.txt --verbose=$(VERBOSE)
	cd /mnt/d/Github/GraphPi/test && /home/yanglaoyuan/.local/lib/python3.8/site-packages/cmake/data/bin/cmake -D TEST_TARGET=asgs_test -D TEST_EXECUTABLE=/mnt/d/Github/GraphPi/bin/asgs_test -D TEST_EXECUTOR= -D TEST_WORKING_DIR=/mnt/d/Github/GraphPi/test -D TEST_EXTRA_ARGS= -D TEST_PROPERTIES= -D TEST_PREFIX= -D TEST_SUFFIX= -D TEST_FILTER= -D NO_PRETTY_TYPES=FALSE -D NO_PRETTY_VALUES=FALSE -D TEST_LIST=asgs_test_TESTS -D CTEST_FILE=/mnt/d/Github/GraphPi/test/asgs_test[1]_tests.cmake -D TEST_DISCOVERY_TIMEOUT=5 -D TEST_XML_OUTPUT_DIR= -P /home/yanglaoyuan/.local/lib/python3.8/site-packages/cmake/data/share/cmake-3.26/Modules/GoogleTestAddTests.cmake

# Rule to build all files generated by this target.
test/CMakeFiles/asgs_test.dir/build: bin/asgs_test
.PHONY : test/CMakeFiles/asgs_test.dir/build

test/CMakeFiles/asgs_test.dir/clean:
	cd /mnt/d/Github/GraphPi/test && $(CMAKE_COMMAND) -P CMakeFiles/asgs_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/asgs_test.dir/clean

test/CMakeFiles/asgs_test.dir/depend:
	cd /mnt/d/Github/GraphPi && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/d/Github/GraphPi /mnt/d/Github/GraphPi/test /mnt/d/Github/GraphPi /mnt/d/Github/GraphPi/test /mnt/d/Github/GraphPi/test/CMakeFiles/asgs_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/asgs_test.dir/depend

