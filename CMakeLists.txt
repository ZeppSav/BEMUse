cmake_minimum_required(VERSION 3.10)

project(BEMUse VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)                                      # Set c++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)                             # Enforce
set(CMAKE_CXX_FLAGS "-std=gnu++11 -O3 -march=native -fopenmp")  # Set compiler optimisation flags

find_package(OpenMP REQUIRED)

#--- Preprocessor Directives
add_compile_definitions(SinglePrec)     # Single precision compile
#add_compile_definitions(DoublePrec)    # Double precision compile

#--- Specify eigen directory
#--- Specify your correct Pfad here!
set(EIGEN_DIR C:\\Development\\eigen\\)

# Specify target executable
add_executable(BEMUse main.cpp)

target_link_libraries(BEMUse PRIVATE OpenMP::OpenMP_CXX)

target_sources(BEMUse
                PRIVATE
                # Headers
                "${PROJECT_SOURCE_DIR}/src/Boundary/FOWT_Platforms.h"
                "${PROJECT_SOURCE_DIR}/src/Geometry/Node.h"

                # Sources
                "${PROJECT_SOURCE_DIR}/src/BEMUser_Console.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Barge.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Boundary_Base.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Ellipsoid.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/GDF_Geo.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Half_Volume_of_Revolution.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Hoizontal_Volume_of_Revolution.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/MAR_Geo.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/STL_Geo.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Ship_Hulls.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Thin_Disc.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Triple_Spar.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Volume_of_Revolution.cpp"
                "${PROJECT_SOURCE_DIR}/src/Boundary/Wing.cpp"
                "${PROJECT_SOURCE_DIR}/src/Geometry/Geo_Elements.cpp"
                "${PROJECT_SOURCE_DIR}/src/Geometry/Panel.cpp"
                "${PROJECT_SOURCE_DIR}/src/Kernels.cpp"
                "${PROJECT_SOURCE_DIR}/src/Solver/Aerodynamic_Solver.cpp"
                "${PROJECT_SOURCE_DIR}/src/Solver/Hydrodynamic_Solver.cpp"
                "${PROJECT_SOURCE_DIR}/src/Solver/Hydrodynamic_Solver_IO.cpp"
                "${PROJECT_SOURCE_DIR}/src/Solver/Solver_Base.cpp"
                "${PROJECT_SOURCE_DIR}/src/Kernels.cpp"
                )

target_include_directories(BEMUse PRIVATE ${EIGEN_DIR})
