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

add_subdirectory(src) # Add the main BEMUse library source directory

# Specify target executable
add_executable(BEMUse src/BEMUser_Console.cpp main.cpp)

target_link_libraries(BEMUse PRIVATE BEMUse_lib)

if(BUILD_GUI)
    add_subdirectory(src_gui)  # Add GUI subdirectory if BUILD_GUI is ON
endif()