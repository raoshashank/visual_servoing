# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /home/shashank/catkin_ws/src/vs_demo_swagat

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default

# Include any dependencies generated for this target.
include CMakeFiles/vs_demo_swagat.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vs_demo_swagat.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vs_demo_swagat.dir/flags.make

CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o: CMakeFiles/vs_demo_swagat.dir/flags.make
CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o: /home/shashank/catkin_ws/src/vs_demo_swagat/src/vsdemo.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o -c /home/shashank/catkin_ws/src/vs_demo_swagat/src/vsdemo.cpp

CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shashank/catkin_ws/src/vs_demo_swagat/src/vsdemo.cpp > CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.i

CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shashank/catkin_ws/src/vs_demo_swagat/src/vsdemo.cpp -o CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.s

CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.requires:
.PHONY : CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.requires

CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.provides: CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.requires
	$(MAKE) -f CMakeFiles/vs_demo_swagat.dir/build.make CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.provides.build
.PHONY : CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.provides

CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.provides.build: CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o

CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o: CMakeFiles/vs_demo_swagat.dir/flags.make
CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o: /home/shashank/catkin_ws/src/vs_demo_swagat/src/ur5model.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o -c /home/shashank/catkin_ws/src/vs_demo_swagat/src/ur5model.cpp

CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shashank/catkin_ws/src/vs_demo_swagat/src/ur5model.cpp > CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.i

CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shashank/catkin_ws/src/vs_demo_swagat/src/ur5model.cpp -o CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.s

CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.requires:
.PHONY : CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.requires

CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.provides: CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.requires
	$(MAKE) -f CMakeFiles/vs_demo_swagat.dir/build.make CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.provides.build
.PHONY : CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.provides

CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.provides.build: CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o

# Object files for target vs_demo_swagat
vs_demo_swagat_OBJECTS = \
"CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o" \
"CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o"

# External object files for target vs_demo_swagat
vs_demo_swagat_EXTERNAL_OBJECTS =

devel/lib/vs_demo_swagat/vs_demo_swagat: CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o
devel/lib/vs_demo_swagat/vs_demo_swagat: CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o
devel/lib/vs_demo_swagat/vs_demo_swagat: CMakeFiles/vs_demo_swagat.dir/build.make
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/libroscpp.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libboost_signals.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/librosconsole.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/librosconsole_log4cxx.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/librosconsole_backend_interface.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/liblog4cxx.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libboost_regex.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/libxmlrpcpp.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/libroscpp_serialization.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/librostime.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /opt/ros/indigo/lib/libcpp_common.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libboost_system.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libboost_thread.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libpthread.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libconsole_bridge.so
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_videostab.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_video.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_superres.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_stitching.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_photo.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_ocl.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_objdetect.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_ml.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_legacy.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_gpu.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_flann.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_core.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_contrib.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /home/shashank/catkin_ws/src/vs_demo_swagat/gnuplot_ci/lib/libgnuplot_ci.a
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_photo.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_legacy.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_video.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_objdetect.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_ml.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_flann.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: /usr/lib/x86_64-linux-gnu/libopencv_core.so.2.4.8
devel/lib/vs_demo_swagat/vs_demo_swagat: CMakeFiles/vs_demo_swagat.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable devel/lib/vs_demo_swagat/vs_demo_swagat"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vs_demo_swagat.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vs_demo_swagat.dir/build: devel/lib/vs_demo_swagat/vs_demo_swagat
.PHONY : CMakeFiles/vs_demo_swagat.dir/build

CMakeFiles/vs_demo_swagat.dir/requires: CMakeFiles/vs_demo_swagat.dir/src/vsdemo.cpp.o.requires
CMakeFiles/vs_demo_swagat.dir/requires: CMakeFiles/vs_demo_swagat.dir/src/ur5model.cpp.o.requires
.PHONY : CMakeFiles/vs_demo_swagat.dir/requires

CMakeFiles/vs_demo_swagat.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vs_demo_swagat.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vs_demo_swagat.dir/clean

CMakeFiles/vs_demo_swagat.dir/depend:
	cd /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shashank/catkin_ws/src/vs_demo_swagat /home/shashank/catkin_ws/src/vs_demo_swagat /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default /home/shashank/catkin_ws/src/build-vs_demo_swagat-Desktop-Default/CMakeFiles/vs_demo_swagat.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vs_demo_swagat.dir/depend

