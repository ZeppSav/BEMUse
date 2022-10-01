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

    -> This geometry object represents a thin disc.
    -> This models a heave plate.

*****************************************************************************/


#ifndef THIN_DISC_H
#define THIN_DISC_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Thin_Disc : public Boundary
{
protected:

    //--- Geo params
    int NR = 0;     // Radial discretisation
    int NA = 0;     // Azimuthal discretisation
    Real R = 0;     // Radius
    Real Depth = 0; // Depth

    //--- Node access
    int Node_ID(int A, int Z);
//    virtual int Aux_Node_ID(int A, int Z)       {}
//    virtual int Ext_Node_ID(int A, int Z)       {}

public:

    //--- Constructor
    Thin_Disc()  {}

    //--- Geometry generation
    void Generate_Nodes();
    void Generate_Elements();
//    virtual void Generate_Aux_Nodes()       {}
//    virtual void Generate_Aux_Elements()    {}
//    virtual void Generate_Ext_Nodes()       {}
//    virtual void Generate_Ext_Elements()    {}

//    void Generate_Symmetry_Nodes(PLANE P);
//    void Generate_Symmetry_Aux_Nodes(PLANE P);

//    void Generate_Symmetry_Elements();
//    void Generate_Symmetry_Aux_Elements();

    //--- Unification of multiple geos
//    virtual void Unify_Geometries()         {}

    //--- Geometry specification
//    virtual void Set_GlobalCoordinateystem(CoordSys *CS)            {Global_CS = CS;}
//    virtual void Set_Coordinateystem(CoordSys *CS)                  {Inertial_CS = CS;}
//    virtual void Set_CoordinateSystems(CoordSys *CSG, CoordSys *CSL) {Global_CS = CSG; Inertial_CS = CSL;}
//    virtual void Set_Origin(Vector3 &O)                             {Origin = O;}

    void Set_Discretisation(std::vector<int> &D)            {NR = D[0]; NA = D[1];}
//    virtual void Set_Auxiliary_Discretisation(std::vector<int> &D)  {}
//    virtual void Set_External_Discretisation(std::vector<int> &D)   {}

    void Set_Dimensions(std::vector<Real> &D)               {R = D[0]; Depth = -D[1];}
//    virtual void Set_Auxiliary_Dimensions(std::vector<Real> &D)     {}
//    virtual void Set_External_Dimensions(std::vector<Real> &D)      {}

//    virtual void Set_Flags(std::vector<bool> &D)                    {Cosine = D[0];}
//    virtual void Set_Nodes(StateVector &D)                          {}

//    //--- Geo Properties
//    Real Volume = 0;
//    Real Surface_Area = 0;
//    Real Max_Dim = 0, Max_Aux_Dim = 0;

//    //--- Kinematic Properties
//    Real Mass = 0;
//    Matrix Mass_Mat = Matrix::Zero(6,6);
//    Vector3 Centre_Gravity = Vector3::Zero();
//    Vector3 Centre_Buoyancy = Vector3::Zero();
//    Matrix3 Radii_Gyration = Matrix3::Zero();

//    //--- Mechanical Properties
//    Matrix Stiffness_Mat = Matrix::Zero(6,6);

//    //--- Geometry calculations
//    virtual void Calculate_GridParams();
//    virtual void Calculate_Volume();
//    virtual void Calculate_COB();
//    virtual void Calculate_SurfaceArea();
//    virtual void Set_Kin_Params()   {}

//    //--- Visualisation Arrays
////    StdVector Vis_Points;
//    std::vector<CMatrix> RadSolArray,DiffSolArray;
//    std::vector<CMatrix> Aux_Vis_Array;

    //--- Destructor
//    ~Boundary();
};

}

#endif // THIN_DISC_H
