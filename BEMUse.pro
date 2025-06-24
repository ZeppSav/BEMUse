#-------------------------------------------------
#
# Project created by QtCreator 2019-03-04T18:40:37
#
#-------------------------------------------------

TARGET = BEMUse

#----------------------------------------------
# Configuration flags
#----------------------------------------------
CONFIG += c++20

#---Optimise flags

QMAKE_CXXFLAGS += -O3               # Optimisations for eigen
QMAKE_CXXFLAGS += -march=native     # Activate all optimisation flags

#----------------------------------------------
# OpenMP support
#----------------------------------------------

QMAKE_CXXFLAGS += -fopenmp  # Compiler support for openMP
LIBS += -fopenmp

#----------------------------------------------
# Specify floating point precision
#----------------------------------------------

DEFINES += SinglePrec
# DEFINES += DoublePrec

#----------------------------------------------
# include path to external source
#----------------------------------------------

INCLUDEPATH += ..\eigen

INCLUDEPATH += src

#------------------------
# Source and header files
#------------------------

SOURCES += main.cpp\
    src/BEMUser_Console.cpp \
    src/Boundary/Barge.cpp \
    src/Boundary/Boundary_Base.cpp \
    src/Boundary/Ellipsoid.cpp \
    src/Boundary/GDF_Geo.cpp \
    src/Boundary/Half_Volume_of_Revolution.cpp \
    src/Boundary/Horizontal_Volume_of_Revolution.cpp \
    src/Boundary/MAR_Geo.cpp \
    src/Boundary/PNL_Geo.cpp \
    src/Boundary/STL_Geo.cpp \
    src/Boundary/Ship_Hulls.cpp \
    src/Boundary/Thin_Disc.cpp \
    src/Boundary/Triple_Spar.cpp \
    src/Boundary/Volume_of_Revolution.cpp \
    src/Boundary/Wing.cpp \
    src/Geometry/Geo_Elements.cpp \
    src/Geometry/Panel.cpp \
    src/Kernels.cpp \
    src/Solver/Aerodynamic_Solver.cpp \
    src/Solver/Hydrodynamic_Solver.cpp \
    src/Solver/Hydrodynamic_Solver_IO.cpp \
    src/Solver/Solver_Base.cpp \

HEADERS  += \
    src/BEMUse_Inquiry.h \
    src/BEMUser_Console.h \
    src/Boundary/Airfoil_Profiles.h \
    src/Boundary/Barge.h \
    src/Boundary/Boundary_Base.h \
    src/Boundary/Ellipsoid.h \
    src/Boundary/FOWT_Platforms.h \
    src/Boundary/GDF_Geo.h \
    src/Boundary/Half_Volume_of_Revolution.h \
    src/Boundary/Horizontal_Volume_of_Revolution.h \
    src/Boundary/MAR_Geo.h \
    src/Boundary/PNL_Geo.h \
    src/Boundary/STL_Geo.h \
    src/Boundary/Ship_Hulls.h \
    src/Boundary/Thin_Disc.h \
    src/Boundary/Triple_Spar.h \
    src/Boundary/Volume_of_Revolution.h \
    src/Boundary/Wing.h \
    src/Geometry/CoordSys.h \
    src/Geometry/Geo_Elements.h \
    src/Geometry/Node.h \
    src/Geometry/Panel.h \
    src/Solver/Solver_Base.h \
    src/Solver/Aerodynamic_Solver.h \
    src/Solver/Hydrodynamic_Solver.h \

#-------------------
# GUI configuration
#-------------------

DEFINES += BEMUse_GUI   # Compile the GUI
QT       += gui opengl widgets openglwidgets

LIBS += -lOpengl32      # include openGL library for visualisation

SOURCES += src_gui/BEMUser_Interface.cpp \
    src_gui/Grid_Options.cpp \
    src_gui/Solver_Setup.cpp \

HEADERS  += src_gui/BEMUser_Interface.h \
    src_gui/Grid_Options.h \
    src_gui/Solver_Setup.h \
    src_gui/Visualise.h \

FORMS    += src_gui/bemuser_interface.ui \
    src_gui/grid_options.ui \
    src_gui/solver_setup.ui
