cmake_minimum_required(VERSION 3.10)

project(BEMUse VERSION 1.0 LANGUAGES CXX)

option(BUILD_GUI "Build GUI version" OFF)  # Option to build GUI version# Optional GUI
option(BUILD_SHARED_LIBS "Build shared library" ON)

IF(NOT DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

set(CMAKE_CXX_STANDARD 17)                     # Set c++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)            # Enforce
set(CMAKE_CXX_FLAGS "-march=native")  # Set compiler optimisation flags

#--- Preprocessor Directives
add_compile_definitions(SinglePrec)     # Single precision compile
#add_compile_definitions(DoublePrec)    # Double precision compile

#--- RPATH tweaks for relocatable installation
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)  # Also add external library paths to RPATH
set(CMAKE_SKIP_BUILD_RPATH FALSE)            # Do not skip build RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)    # Use install RPATH during build
# Set RPATH to point to the lib directory relative to the executable
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin") # Not tested
    list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
endif()
# On Windows, shared libraries are typically in the same directory as the executable

#--- Add library directory
add_subdirectory(src) # Add the main BEMUse library source directory

# Specify target executable
add_executable(BEMUse src/BEMUser_Console.cpp main.cpp)

target_link_libraries(BEMUse PRIVATE BEMUse_lib)

install(TARGETS BEMUse)

if(BUILD_GUI)
    add_subdirectory(src_gui)  # Add GUI subdirectory if BUILD_GUI is ON
endif()