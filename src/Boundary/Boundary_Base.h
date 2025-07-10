/****************************************************************************
    BEMUse Boundary
    Copyright (C) 2022 Joseph Saverin j.saverin@tu-berlin.de

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Information on file:

    -> This is the header which contains the BEMUse Boundary type

*****************************************************************************/

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../BEMUse_IO.h"
#include "../Geometry/Geo_Elements.h"

namespace BEMUse
{

// This is the base class for a boundary. It stores the nodes and elements required for the analysis.
// Generation of panel elements depends upon the type of solver being applied, so this is avoided in this class.

//--- Discretisation functions

inline Real PiCosFac(int i, int itot)      {return 0.5*(1.0-cos(i*PI/(itot-1)));}
inline Real HalfPiCosFac(int i, int itot)  {return -cos(0.5*PI+0.5*i*PI/(itot-1));}
inline Real HalfPiCosFac2(int i, int itot)  {return 1.0 - cos(0.5*i*PI/(itot-1));}

inline Real PiCosFacShift(int i, int itot) {return 0.5*(1.0-cos((i+0.5)*PI/(itot-1)));}
inline Real HalfPiCosFacShift(int i, int itot) {return -cos(0.5*PI+0.5*(i+0.5)*PI/(itot-1));}
inline Real HalfPiCosFac2Shift(int i, int itot)  {return 1.0 - cos(0.5*(i+0.5)*PI/(itot-1));}

struct Parameter : private std::tuple<std::string, int, Real, bool>
{
    using Base = std::tuple<std::string, int, Real, bool>;
    Parameter(const std::string& name, int val)   : Base(name, val, Real{}, bool{}) {}
    Parameter(const std::string& name, Real val)  : Base(name, int{}, val, bool{})  {}
    Parameter(const std::string& name, bool val)  : Base(name, int{}, Real{}, val)  {}
    bool myNameis(const std::string &a)    {return (a == std::get<0>(*this));}
    template <typename T> T Get_Param();    // Parameter retrieval: Templ. fns defined in Boundary_Base.cpp
};

enum PLANE  {XPlane, YPlane, ZPlane};

class Boundary
{

protected:

    //--- Parameters for geometry
    std::vector<Parameter> Parameters;

    //--- Nodes
    std::vector<SP_Node>    Nodes,          Symm_Nodes;
    std::vector<SP_Node>    Aux_Nodes,      Symm_Aux_Nodes;
    std::vector<SP_Node>    FreeSurface_Nodes;
    std::vector<SP_Node>    Wall_Nodes;
    std::vector<SP_Node>    Floor_Nodes;

    //--- Elements
    std::vector<SP_Geo>     Elements,       Symm_Elements;
    std::vector<SP_Geo>     Aux_Elements,   Symm_Aux_Elements;
    std::vector<SP_Geo>     FreeSurface_Elements;
    std::vector<SP_Geo>     Wall_Elements;
    std::vector<SP_Geo>     Floor_Elements;

    bool VisWall = true;
    bool VisFloor = true;

    //--- Geo flags
    bool TriPanels = false;
    bool Cosine = true;
    bool Panels_Outward = false;

    //--- Node access
    virtual int Node_ID(int A, int Z)           {}
    virtual int Aux_Node_ID(int A, int Z)       {}
    virtual int Ext_Node_ID(int A, int Z)       {}
    virtual int Wall_ID(int A, int Z)           {}

    CoordSys *Global_CS = nullptr;
    CoordSys *Inertial_CS = nullptr;

public:

    //--- Constructor
    Boundary()  {
        // if (!Global_CS)    Global_CS = new CoordSys();
        // if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);
    }

    //--- Geo setup
    virtual void Setup();
    virtual void Set_Parameters(std::vector<Parameter> &Params)       {}

    //--- Import/Export
    virtual bool Read_Input_File(std::string &FilePath) {}
    void Export_Geometry_GDF(std::string &FilePath);
    void Export_Geometry_STL(std::string &FilePath);
    void Export_Geometry_MAR(std::string &FilePath);
    std::string ImportFile = "", FileInfo = "";

    //--- Geometry generation
    virtual void Generate_Nodes()           {}
    virtual void Generate_Elements()        {}
    virtual void Generate_Aux_Nodes()       {}
    virtual void Generate_Aux_Elements()    {}
    virtual void Generate_FreeSurface_Nodes()       {}
    virtual void Generate_FreeSurface_Elements()    {}
    virtual void Generate_Wall_Nodes()      {}
    virtual void Generate_Wall_Elements()   {}
    virtual void Generate_Floor_Nodes()     {}
    virtual void Generate_Floor_Elements()  {}

    void Generate_Symmetry_Nodes(PLANE P);
    void Generate_Symmetry_Aux_Nodes(PLANE P);

    void Generate_Symmetry_Elements();
    void Generate_Symmetry_Aux_Elements();

    //--- Unification of multiple geos
    virtual void Unify_Geometries()         {}

    //--- Geometry specification
    //    virtual void Set_GlobalCoordinateystem(CoordSys *CS)            {Global_CS = CS;}
    //    virtual void Set_Coordinateystem(CoordSys *CS)                  {Inertial_CS = CS;}
    virtual void Set_CoordinateSystems(CoordSys *CSG, CoordSys *CSL) {Global_CS = CSG; Inertial_CS = CSL;}

    //--- Geo Properties
    Real Volume = 0;
    Real Surface_Area = 0;
    Real Max_Dim = 0, Max_Aux_Dim = 0;

    //--- Kinematic Properties
    Real Mass = 0;
    Matrix Mass_Mat = Matrix::Zero(6,6);
    Vector3 Centre_Gravity = Vector3::Zero();
    Vector3 Centre_Buoyancy = Vector3::Zero();
    Matrix3 Radii_Gyration = Matrix3::Zero();

    //--- Mechanical Properties
    Matrix Stiffness_Mat = Matrix::Zero(6,6);

    //--- Geometry calculations
    virtual void Calculate_GridParams();
    virtual void Calculate_Volume();
    virtual void Calculate_COB();
    virtual void Calculate_SurfaceArea();
    virtual void Set_Kin_Params()   {}

    //--- Retrieve grid stats
    int N_Nodes()           {return Nodes.size();}
    int N_Aux_Nodes()       {return Aux_Nodes.size();}
    int N_Elements()        {return Elements.size();}
    int N_Aux_Elements()    {return Aux_Elements.size();}
    int N_WallRDim()        {return 0;}
    int N_WallAziDim()      {return 0;}

    //--- Retrieve elements for visualisation
    void    Get_Nodes(std::vector<SP_Node> &L)              {StdAppend(L,Nodes);}
    void    Get_Elements(std::vector<SP_Geo> &L)            {StdAppend(L,Elements); }
    void    Get_Symm_Elements(std::vector<SP_Geo> &L)       {StdAppend(L,Symm_Elements);}
    void    Get_Aux_Nodes(std::vector<SP_Node> &L)          {StdAppend(L,Aux_Nodes);}
    void    Get_Aux_Elements(std::vector<SP_Geo> &L)        {StdAppend(L,Aux_Elements);}
    void    Get_Aux_Symm_Elements(std::vector<SP_Geo> &L)   {StdAppend(L,Symm_Aux_Elements);}
    void    Get_FreeSurface_Nodes(std::vector<SP_Node> &L)  {StdAppend(L,FreeSurface_Nodes);}
    void    Get_FreeSurface_Elements(std::vector<SP_Geo> &L){StdAppend(L,FreeSurface_Elements);}
    void    Get_Floor_Elements(std::vector<SP_Geo> &L)      {StdAppend(L,Floor_Elements);}
    void    Get_Wall_Elements(std::vector<SP_Geo> &L)       {StdAppend(L,Wall_Elements);}
    void    Get_Wall_Nodes(std::vector<SP_Node> &L)         {StdAppend(L,Wall_Nodes);}
    void    Get_Floor_Nodes(std::vector<SP_Node> &L)        {StdAppend(L,Floor_Nodes);}
    void    Get_SurfaceWall_Elements(std::vector<SP_Geo> &L)     {
        StdAppend(L,Elements);
        if (VisWall)    StdAppend(L,Wall_Elements);
        if (VisFloor)   StdAppend(L,Floor_Elements);
    }

    CoordSys *Get_Inertial_CS()                             {return Inertial_CS;}

    //--- Visualisation Arrays
    //    StdVector Vis_Points;
    std::vector<CMatrix> RadSolArray;
    std::vector<CMatrix> DiffSolArray;
    std::vector<CMatrix> FS_Rad_Array;
    std::vector<CMatrix> FS_Scat_Array;
    void Toggle_Wall_Vis(bool f)    {VisWall = f;}
    void Toggle_Floor_Vis(bool f)   {VisFloor = f;}

    //--- Retrieve parameters
    template <class T>
    inline bool Retrieve_Parameter(const std::string &a, T &val)
    {
        for (Parameter P : Parameters) {
            if (P.myNameis(a)) {
                val = P.Get_Param<T>();
                return true;
            }
        }
        return false;
    }

    //--- Destructor
    virtual ~Boundary();
};

}

#endif // BOUNDARY_H
