# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build

# Include any dependencies generated for this target.
include src/swig/CMakeFiles/_nlopt.dir/depend.make

# Include the progress variables for this target.
include src/swig/CMakeFiles/_nlopt.dir/progress.make

# Include the compile flags for this target's objects.
include src/swig/CMakeFiles/_nlopt.dir/flags.make

src/swig/CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.o: src/swig/CMakeFiles/_nlopt.dir/flags.make
src/swig/CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.o: src/swig/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/swig/CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.o"
	cd /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.o -c /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx

src/swig/CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.i"
	cd /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx > CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.i

src/swig/CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.s"
	cd /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx -o CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.s

# Object files for target _nlopt
_nlopt_OBJECTS = \
"CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.o"

# External object files for target _nlopt
_nlopt_EXTERNAL_OBJECTS =

src/swig/_nlopt.so: src/swig/CMakeFiles/_nlopt.dir/CMakeFiles/_nlopt.dir/nloptPYTHON_wrap.cxx.o
src/swig/_nlopt.so: src/swig/CMakeFiles/_nlopt.dir/build.make
src/swig/_nlopt.so: libnlopt.a
src/swig/_nlopt.so: /usr/lib/libpython3.8.so
src/swig/_nlopt.so: src/swig/CMakeFiles/_nlopt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module _nlopt.so"
	cd /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_nlopt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/swig/CMakeFiles/_nlopt.dir/build: src/swig/_nlopt.so

.PHONY : src/swig/CMakeFiles/_nlopt.dir/build

src/swig/CMakeFiles/_nlopt.dir/clean:
	cd /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig && $(CMAKE_COMMAND) -P CMakeFiles/_nlopt.dir/cmake_clean.cmake
.PHONY : src/swig/CMakeFiles/_nlopt.dir/clean

src/swig/CMakeFiles/_nlopt.dir/depend:
	cd /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0 /home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0/src/swig /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig /media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out/build/src/swig/CMakeFiles/_nlopt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/swig/CMakeFiles/_nlopt.dir/depend
