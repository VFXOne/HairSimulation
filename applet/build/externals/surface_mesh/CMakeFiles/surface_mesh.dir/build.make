# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = "/home/lucas/Documents/Semester project/proj_dyn"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/lucas/Documents/Semester project/proj_dyn/build"

# Include any dependencies generated for this target.
include externals/surface_mesh/CMakeFiles/surface_mesh.dir/depend.make

# Include the progress variables for this target.
include externals/surface_mesh/CMakeFiles/surface_mesh.dir/progress.make

# Include the compile flags for this target's objects.
include externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o: externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make
externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o: ../externals/surface_mesh/surface_mesh/IO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o -c "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO.cpp"

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.i"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO.cpp" > CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.i

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.s"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO.cpp" -o CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.s

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.requires:

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.provides: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.requires
	$(MAKE) -f externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.provides.build
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.provides

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.provides.build: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o


externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o: externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make
externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o: ../externals/surface_mesh/surface_mesh/IO_obj.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o -c "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_obj.cpp"

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.i"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_obj.cpp" > CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.i

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.s"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_obj.cpp" -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.s

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.requires:

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.provides: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.requires
	$(MAKE) -f externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.provides.build
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.provides

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.provides.build: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o


externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o: externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make
externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o: ../externals/surface_mesh/surface_mesh/IO_off.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o -c "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_off.cpp"

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.i"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_off.cpp" > CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.i

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.s"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_off.cpp" -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.s

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.requires:

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.provides: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.requires
	$(MAKE) -f externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.provides.build
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.provides

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.provides.build: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o


externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o: externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make
externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o: ../externals/surface_mesh/surface_mesh/IO_poly.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o -c "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_poly.cpp"

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.i"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_poly.cpp" > CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.i

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.s"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_poly.cpp" -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.s

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.requires:

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.provides: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.requires
	$(MAKE) -f externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.provides.build
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.provides

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.provides.build: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o


externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o: externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make
externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o: ../externals/surface_mesh/surface_mesh/IO_stl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o -c "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_stl.cpp"

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.i"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_stl.cpp" > CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.i

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.s"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/IO_stl.cpp" -o CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.s

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.requires:

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.provides: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.requires
	$(MAKE) -f externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.provides.build
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.provides

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.provides.build: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o


externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o: externals/surface_mesh/CMakeFiles/surface_mesh.dir/flags.make
externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o: ../externals/surface_mesh/surface_mesh/Surface_mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o -c "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/Surface_mesh.cpp"

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.i"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/Surface_mesh.cpp" > CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.i

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.s"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh/surface_mesh/Surface_mesh.cpp" -o CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.s

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.requires:

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.provides: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.requires
	$(MAKE) -f externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.provides.build
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.provides

externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.provides.build: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o


# Object files for target surface_mesh
surface_mesh_OBJECTS = \
"CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o" \
"CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o" \
"CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o" \
"CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o" \
"CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o" \
"CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o"

# External object files for target surface_mesh
surface_mesh_EXTERNAL_OBJECTS =

externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/build.make
externals/surface_mesh/libsurface_mesh.so.1.0: externals/surface_mesh/CMakeFiles/surface_mesh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/lucas/Documents/Semester project/proj_dyn/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX shared library libsurface_mesh.so"
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/surface_mesh.dir/link.txt --verbose=$(VERBOSE)
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && $(CMAKE_COMMAND) -E cmake_symlink_library libsurface_mesh.so.1.0 libsurface_mesh.so.1.0 libsurface_mesh.so

externals/surface_mesh/libsurface_mesh.so: externals/surface_mesh/libsurface_mesh.so.1.0
	@$(CMAKE_COMMAND) -E touch_nocreate externals/surface_mesh/libsurface_mesh.so

# Rule to build all files generated by this target.
externals/surface_mesh/CMakeFiles/surface_mesh.dir/build: externals/surface_mesh/libsurface_mesh.so

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/build

externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO.cpp.o.requires
externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_obj.cpp.o.requires
externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_off.cpp.o.requires
externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_poly.cpp.o.requires
externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/IO_stl.cpp.o.requires
externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires: externals/surface_mesh/CMakeFiles/surface_mesh.dir/surface_mesh/Surface_mesh.cpp.o.requires

.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/requires

externals/surface_mesh/CMakeFiles/surface_mesh.dir/clean:
	cd "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" && $(CMAKE_COMMAND) -P CMakeFiles/surface_mesh.dir/cmake_clean.cmake
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/clean

externals/surface_mesh/CMakeFiles/surface_mesh.dir/depend:
	cd "/home/lucas/Documents/Semester project/proj_dyn/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/lucas/Documents/Semester project/proj_dyn" "/home/lucas/Documents/Semester project/proj_dyn/externals/surface_mesh" "/home/lucas/Documents/Semester project/proj_dyn/build" "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh" "/home/lucas/Documents/Semester project/proj_dyn/build/externals/surface_mesh/CMakeFiles/surface_mesh.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : externals/surface_mesh/CMakeFiles/surface_mesh.dir/depend

