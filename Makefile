# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony/CMakeFiles /home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named UserSubLib

# Build rule for target.
UserSubLib: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserSubLib
.PHONY : UserSubLib

# fast build rule for target.
UserSubLib/fast:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/build
.PHONY : UserSubLib/fast

UELMAT_Assembly.o: UELMAT_Assembly.cpp.o

.PHONY : UELMAT_Assembly.o

# target to build an object file
UELMAT_Assembly.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Assembly.cpp.o
.PHONY : UELMAT_Assembly.cpp.o

UELMAT_Assembly.i: UELMAT_Assembly.cpp.i

.PHONY : UELMAT_Assembly.i

# target to preprocess a source file
UELMAT_Assembly.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Assembly.cpp.i
.PHONY : UELMAT_Assembly.cpp.i

UELMAT_Assembly.s: UELMAT_Assembly.cpp.s

.PHONY : UELMAT_Assembly.s

# target to generate assembly for a file
UELMAT_Assembly.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Assembly.cpp.s
.PHONY : UELMAT_Assembly.cpp.s

UELMAT_Integration.o: UELMAT_Integration.cpp.o

.PHONY : UELMAT_Integration.o

# target to build an object file
UELMAT_Integration.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Integration.cpp.o
.PHONY : UELMAT_Integration.cpp.o

UELMAT_Integration.i: UELMAT_Integration.cpp.i

.PHONY : UELMAT_Integration.i

# target to preprocess a source file
UELMAT_Integration.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Integration.cpp.i
.PHONY : UELMAT_Integration.cpp.i

UELMAT_Integration.s: UELMAT_Integration.cpp.s

.PHONY : UELMAT_Integration.s

# target to generate assembly for a file
UELMAT_Integration.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Integration.cpp.s
.PHONY : UELMAT_Integration.cpp.s

UELMAT_Matrices.o: UELMAT_Matrices.cpp.o

.PHONY : UELMAT_Matrices.o

# target to build an object file
UELMAT_Matrices.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Matrices.cpp.o
.PHONY : UELMAT_Matrices.cpp.o

UELMAT_Matrices.i: UELMAT_Matrices.cpp.i

.PHONY : UELMAT_Matrices.i

# target to preprocess a source file
UELMAT_Matrices.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Matrices.cpp.i
.PHONY : UELMAT_Matrices.cpp.i

UELMAT_Matrices.s: UELMAT_Matrices.cpp.s

.PHONY : UELMAT_Matrices.s

# target to generate assembly for a file
UELMAT_Matrices.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Matrices.cpp.s
.PHONY : UELMAT_Matrices.cpp.s

UELMAT_Order.o: UELMAT_Order.cpp.o

.PHONY : UELMAT_Order.o

# target to build an object file
UELMAT_Order.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Order.cpp.o
.PHONY : UELMAT_Order.cpp.o

UELMAT_Order.i: UELMAT_Order.cpp.i

.PHONY : UELMAT_Order.i

# target to preprocess a source file
UELMAT_Order.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Order.cpp.i
.PHONY : UELMAT_Order.cpp.i

UELMAT_Order.s: UELMAT_Order.cpp.s

.PHONY : UELMAT_Order.s

# target to generate assembly for a file
UELMAT_Order.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_Order.cpp.s
.PHONY : UELMAT_Order.cpp.s

UELMAT_RTLM.o: UELMAT_RTLM.cpp.o

.PHONY : UELMAT_RTLM.o

# target to build an object file
UELMAT_RTLM.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_RTLM.cpp.o
.PHONY : UELMAT_RTLM.cpp.o

UELMAT_RTLM.i: UELMAT_RTLM.cpp.i

.PHONY : UELMAT_RTLM.i

# target to preprocess a source file
UELMAT_RTLM.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_RTLM.cpp.i
.PHONY : UELMAT_RTLM.cpp.i

UELMAT_RTLM.s: UELMAT_RTLM.cpp.s

.PHONY : UELMAT_RTLM.s

# target to generate assembly for a file
UELMAT_RTLM.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_RTLM.cpp.s
.PHONY : UELMAT_RTLM.cpp.s

UELMAT_ShapeFunctions.o: UELMAT_ShapeFunctions.cpp.o

.PHONY : UELMAT_ShapeFunctions.o

# target to build an object file
UELMAT_ShapeFunctions.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_ShapeFunctions.cpp.o
.PHONY : UELMAT_ShapeFunctions.cpp.o

UELMAT_ShapeFunctions.i: UELMAT_ShapeFunctions.cpp.i

.PHONY : UELMAT_ShapeFunctions.i

# target to preprocess a source file
UELMAT_ShapeFunctions.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_ShapeFunctions.cpp.i
.PHONY : UELMAT_ShapeFunctions.cpp.i

UELMAT_ShapeFunctions.s: UELMAT_ShapeFunctions.cpp.s

.PHONY : UELMAT_ShapeFunctions.s

# target to generate assembly for a file
UELMAT_ShapeFunctions.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/UELMAT_ShapeFunctions.cpp.s
.PHONY : UELMAT_ShapeFunctions.cpp.s

USUB_UtilityRoutines.o: USUB_UtilityRoutines.cpp.o

.PHONY : USUB_UtilityRoutines.o

# target to build an object file
USUB_UtilityRoutines.cpp.o:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/USUB_UtilityRoutines.cpp.o
.PHONY : USUB_UtilityRoutines.cpp.o

USUB_UtilityRoutines.i: USUB_UtilityRoutines.cpp.i

.PHONY : USUB_UtilityRoutines.i

# target to preprocess a source file
USUB_UtilityRoutines.cpp.i:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/USUB_UtilityRoutines.cpp.i
.PHONY : USUB_UtilityRoutines.cpp.i

USUB_UtilityRoutines.s: USUB_UtilityRoutines.cpp.s

.PHONY : USUB_UtilityRoutines.s

# target to generate assembly for a file
USUB_UtilityRoutines.cpp.s:
	$(MAKE) -f CMakeFiles/UserSubLib.dir/build.make CMakeFiles/UserSubLib.dir/USUB_UtilityRoutines.cpp.s
.PHONY : USUB_UtilityRoutines.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... install/strip"
	@echo "... edit_cache"
	@echo "... UserSubLib"
	@echo "... rebuild_cache"
	@echo "... list_install_components"
	@echo "... install/local"
	@echo "... install"
	@echo "... UELMAT_Assembly.o"
	@echo "... UELMAT_Assembly.i"
	@echo "... UELMAT_Assembly.s"
	@echo "... UELMAT_Integration.o"
	@echo "... UELMAT_Integration.i"
	@echo "... UELMAT_Integration.s"
	@echo "... UELMAT_Matrices.o"
	@echo "... UELMAT_Matrices.i"
	@echo "... UELMAT_Matrices.s"
	@echo "... UELMAT_Order.o"
	@echo "... UELMAT_Order.i"
	@echo "... UELMAT_Order.s"
	@echo "... UELMAT_RTLM.o"
	@echo "... UELMAT_RTLM.i"
	@echo "... UELMAT_RTLM.s"
	@echo "... UELMAT_ShapeFunctions.o"
	@echo "... UELMAT_ShapeFunctions.i"
	@echo "... UELMAT_ShapeFunctions.s"
	@echo "... USUB_UtilityRoutines.o"
	@echo "... USUB_UtilityRoutines.i"
	@echo "... USUB_UtilityRoutines.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

