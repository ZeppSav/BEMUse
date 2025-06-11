# Dependencies
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

add_library(BEMUse_lib)  # Create a static library for BEMUse
set_target_properties(BEMUse_lib PROPERTIES OUTPUT_NAME BEMUse)

# Sources
target_sources(BEMUse_lib PRIVATE
    Boundary/Barge.cpp
    Boundary/Boundary_Base.cpp
    Boundary/Ellipsoid.cpp
    Boundary/GDF_Geo.cpp
    Boundary/Half_Volume_of_Revolution.cpp
    Boundary/Hoizontal_Volume_of_Revolution.cpp
    Boundary/MAR_Geo.cpp
    Boundary/PNL_Geo.cpp
    Boundary/STL_Geo.cpp
    Boundary/Ship_Hulls.cpp
    Boundary/Thin_Disc.cpp
    Boundary/Triple_Spar.cpp
    Boundary/Volume_of_Revolution.cpp
    Boundary/Wing.cpp
    Geometry/Geo_Elements.cpp
    Geometry/Panel.cpp
    Kernels.cpp
    Solver/Aerodynamic_Solver.cpp
    Solver/Hydrodynamic_Solver.cpp
    Solver/Hydrodynamic_Solver_IO.cpp
    Solver/Solver_Base.cpp
)

# Headers
target_sources(BEMUse_lib PUBLIC
    FILE_SET HEADERS
    FILES
    Boundary/FOWT_Platforms.h
    Geometry/Node.h
    Boundary/Ellipsoid.h
    Boundary/Volume_of_Revolution.h
    Boundary/Horizontal_Volume_of_Revolution.h
    Boundary/Triple_Spar.h
    Boundary/Ship_Hulls.h
    Boundary/STL_Geo.h
    Boundary/GDF_Geo.h
    Boundary/MAR_Geo.h
    Boundary/PNL_Geo.h
    Boundary/Barge.h
    Boundary/Thin_Disc.h
    Boundary/Wing.h
)

target_link_libraries(BEMUse_lib PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS BEMUse_lib FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/BEMUse)