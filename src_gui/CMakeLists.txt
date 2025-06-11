find_package(Qt6 REQUIRED COMPONENTS Gui Widgets OpenGL OpenGLWidgets)
find_package(OpenGL REQUIRED)

add_executable(BEMUse_GUI)

target_sources(BEMUse_GUI PRIVATE
    ../main.cpp
    BEMUser_Interface.cpp
    Grid_Options.cpp
    Solver_Setup.cpp
)

target_sources(BEMUse_GUI PRIVATE
    FILE_SET HEADERS
    FILES
    BEMUser_Interface.h
    Grid_Options.h
    Solver_Setup.h
)

set_target_properties(BEMUse_GUI PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

target_compile_definitions(BEMUse_GUI PRIVATE BEMUse_GUI)

set_target_properties(BEMUse_GUI PROPERTIES AUTOUIC ON AUTOMOC ON)

target_link_libraries(BEMUse_GUI PRIVATE BEMUse_lib Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets OpenGL::GLU)