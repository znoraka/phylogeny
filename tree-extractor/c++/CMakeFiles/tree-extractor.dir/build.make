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
CMAKE_SOURCE_DIR = /home/noe/Documents/dev/phylogeny/tree-extractor/c++

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/noe/Documents/dev/phylogeny/tree-extractor/c++

# Include any dependencies generated for this target.
include CMakeFiles/tree-extractor.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tree-extractor.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tree-extractor.dir/flags.make

CMakeFiles/tree-extractor.dir/src/main.cpp.o: CMakeFiles/tree-extractor.dir/flags.make
CMakeFiles/tree-extractor.dir/src/main.cpp.o: src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/noe/Documents/dev/phylogeny/tree-extractor/c++/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tree-extractor.dir/src/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tree-extractor.dir/src/main.cpp.o -c /home/noe/Documents/dev/phylogeny/tree-extractor/c++/src/main.cpp

CMakeFiles/tree-extractor.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tree-extractor.dir/src/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/noe/Documents/dev/phylogeny/tree-extractor/c++/src/main.cpp > CMakeFiles/tree-extractor.dir/src/main.cpp.i

CMakeFiles/tree-extractor.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tree-extractor.dir/src/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/noe/Documents/dev/phylogeny/tree-extractor/c++/src/main.cpp -o CMakeFiles/tree-extractor.dir/src/main.cpp.s

CMakeFiles/tree-extractor.dir/src/main.cpp.o.requires:

.PHONY : CMakeFiles/tree-extractor.dir/src/main.cpp.o.requires

CMakeFiles/tree-extractor.dir/src/main.cpp.o.provides: CMakeFiles/tree-extractor.dir/src/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/tree-extractor.dir/build.make CMakeFiles/tree-extractor.dir/src/main.cpp.o.provides.build
.PHONY : CMakeFiles/tree-extractor.dir/src/main.cpp.o.provides

CMakeFiles/tree-extractor.dir/src/main.cpp.o.provides.build: CMakeFiles/tree-extractor.dir/src/main.cpp.o


# Object files for target tree-extractor
tree__extractor_OBJECTS = \
"CMakeFiles/tree-extractor.dir/src/main.cpp.o"

# External object files for target tree-extractor
tree__extractor_EXTERNAL_OBJECTS =

bin/tree-extractor: CMakeFiles/tree-extractor.dir/src/main.cpp.o
bin/tree-extractor: CMakeFiles/tree-extractor.dir/build.make
bin/tree-extractor: /usr/lib/x86_64-linux-gnu/libMagick++-6.Q16.so
bin/tree-extractor: /usr/lib/x86_64-linux-gnu/libMagickCore-6.Q16.so
bin/tree-extractor: /usr/lib/x86_64-linux-gnu/libMagickWand-6.Q16.so
bin/tree-extractor: /usr/lib/x86_64-linux-gnu/libjpeg.so
bin/tree-extractor: CMakeFiles/tree-extractor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/noe/Documents/dev/phylogeny/tree-extractor/c++/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/tree-extractor"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tree-extractor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tree-extractor.dir/build: bin/tree-extractor

.PHONY : CMakeFiles/tree-extractor.dir/build

CMakeFiles/tree-extractor.dir/requires: CMakeFiles/tree-extractor.dir/src/main.cpp.o.requires

.PHONY : CMakeFiles/tree-extractor.dir/requires

CMakeFiles/tree-extractor.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tree-extractor.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tree-extractor.dir/clean

CMakeFiles/tree-extractor.dir/depend:
	cd /home/noe/Documents/dev/phylogeny/tree-extractor/c++ && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/noe/Documents/dev/phylogeny/tree-extractor/c++ /home/noe/Documents/dev/phylogeny/tree-extractor/c++ /home/noe/Documents/dev/phylogeny/tree-extractor/c++ /home/noe/Documents/dev/phylogeny/tree-extractor/c++ /home/noe/Documents/dev/phylogeny/tree-extractor/c++/CMakeFiles/tree-extractor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tree-extractor.dir/depend

