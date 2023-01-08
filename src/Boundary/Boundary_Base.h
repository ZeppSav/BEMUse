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

enum PLANE  {XPlane, YPlane, ZPlane};

class Boundary
{

protected:

    //--- Nodes
    std::vector<SP_Node>    Nodes,          Symm_Nodes;
    std::vector<SP_Node>    Aux_Nodes,      Symm_Aux_Nodes,     Ext_Nodes;

    //--- Elements
    std::vector<SP_Geo>     Elements,       Symm_Elements;
    std::vector<SP_Geo>     Aux_Elements,   Symm_Aux_Elements,  Ext_Elements;

    //--- Geo flags
    bool Cosine = true;

    //--- Node access
    virtual int Node_ID(int A, int Z)           {}
    virtual int Aux_Node_ID(int A, int Z)       {}
    virtual int Ext_Node_ID(int A, int Z)       {}

    CoordSys *Global_CS = NULL;
    CoordSys *Inertial_CS = NULL;
//    Vector3 Origin;

public:

    //--- Constructor
    Boundary()  {}

    //--- Geo setup
    virtual void Setup();

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
    virtual void Generate_Ext_Nodes()       {}
    virtual void Generate_Ext_Elements()    {}

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
//    virtual void Set_Origin(Vector3 &O)                             {Origin = O;}

    virtual void Set_Discretisation(std::vector<int> &D)            {}
    virtual void Set_Auxiliary_Discretisation(std::vector<int> &D)  {}
    virtual void Set_External_Discretisation(std::vector<int> &D)   {}

    virtual void Set_Dimensions(std::vector<Real> &D)               {}
    virtual void Set_Auxiliary_Dimensions(std::vector<Real> &D)     {}
    virtual void Set_External_Dimensions(std::vector<Real> &D)      {}

    virtual void Set_Flags(std::vector<bool> &D)                    {Cosine = D[0];}
//    virtual void Set_Nodes(StateVector &D)                          {}

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

    //--- Retrieve elements for visualisation
    void    Get_Nodes(std::vector<SP_Node> &L)              {StdAppend(L,Nodes);}
    void    Get_Elements(std::vector<SP_Geo> &L)            {StdAppend(L,Elements);}
    void    Get_Symm_Elements(std::vector<SP_Geo> &L)       {StdAppend(L,Symm_Elements);}
    void    Get_Aux_Nodes(std::vector<SP_Node> &L)          {StdAppend(L,Aux_Nodes);}
    void    Get_Aux_Elements(std::vector<SP_Geo> &L)        {StdAppend(L,Aux_Elements);}
    void    Get_Aux_Symm_Elements(std::vector<SP_Geo> &L)   {StdAppend(L,Symm_Aux_Elements);}
    void    Get_Ext_Nodes(std::vector<SP_Node> &L)          {StdAppend(L,Ext_Nodes);}
    void    Get_Ext_Elements(std::vector<SP_Geo> &L)        {StdAppend(L,Ext_Elements);}
    CoordSys *Get_Inertial_CS()                             {return Inertial_CS;}

    //--- Visualisation Arrays
//    StdVector Vis_Points;
    std::vector<CMatrix> RadSolArray;
    std::vector<CMatrix> DiffSolArray;
    std::vector<CMatrix> FS_Rad_Array;
    std::vector<CMatrix> FS_Scat_Array;

    //--- Destructor
    ~Boundary();
};

}

#endif // BOUNDARY_H
