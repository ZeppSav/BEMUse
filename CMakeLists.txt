cmake_minimum_required(VERSION 3.10)

project(BEMUse VERSION 1.0 LANGUAGES CXX)

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
