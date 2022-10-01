/****************************************************************************
    BEMUse Volume of revolution
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

    -> This generates a volume of revolution. The variable "Open" specifies whether the
    -> Geometry is a closed volume of revolution or not.

*****************************************************************************/

#ifndef VOLUME_OF_REVOLUTION_H
#define VOLUME_OF_REVOLUTION_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Volume_of_Revolution : public Boundary
{

protected:

    //--- Geo variables
    Real R;

    //--- Perimeter variables
//    Vector Origin = Vector3(0,0,0);
    Real RFS = 0;
    int NA, NZ;
    int NRFS=0, NRES=0;
    std::vector<Vector3> Perimeter;

public:

    //--- Constructor
    Volume_of_Revolution()   {}

    //--- Geometry specification
//    void Set_Origin(Vector3 &O)                     {Origin = O;}
//    void Set_Flags(std::vector<bool> &D)            {}
//    void Set_Nodes(StateVector &D)                  {for (Vector V : D) Perimeter.push_back(Vector3(V(0),V(1),V(2)));}

    //--- Geometry functions
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes();
    void Generate_Aux_Elements();
    void Generate_Ext_Nodes();
    void Generate_Ext_Elements();

    int Node_ID(int A, int Z);
    int Aux_Node_ID(int A, int Z)       {}
    int Ext_Node_ID(int A, int Z);
};

class Half_Cylinder : public Volume_of_Revolution
{

protected:

    //--- Geo variables
    int NR, NV;
    Real Z;

public:

    //--- Constructor
    Half_Cylinder()   {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NA = D[0]; NR = D[1]; NV = D[2];}
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {NRFS = D[0];}
    void Set_External_Discretisation(std::vector<int> &D)   {NRES = D[0];}

    void Set_Dimensions(std::vector<Real> &D)               {R = D[0]; Z = -D[1];}
    void Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0];}

    //--- Geometry functions
    void Generate_Nodes();

};

class Tapered_SparBuoy : public Volume_of_Revolution
{

protected:

    // Perimeter variables
    int NR, NV1, NV2, NV3, NV4;
    Real Z, RB, RT, H1, H2, H3;

public:

    //--- Constructor
    Tapered_SparBuoy()   {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NA = D[0]; NR = D[1]; NV1 = D[2]; NV2 = D[3]; NV3 = D[4]; }
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {NRFS = D[0]; }
    void Set_External_Discretisation(std::vector<int> &D)   {NRES = D[0];}

    void Set_Dimensions(std::vector<Real> &D)               {RB = D[0];  RT = D[1]; Z = -D[2];  H1 = D[3];  H2 = D[4]; }
    void Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0]; R = RT;}

    //--- Geometry functions
    void Generate_Nodes();
};

class Spar_Leg : public Volume_of_Revolution
{

protected:

    // Perimeter variables
    int N1, N2, N3, N4;
    Real Z, RB, RT, H1;

public:

    //--- Constructor
    Spar_Leg()   {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NA = D[0]; N1 = D[1]; N2 = D[2]; N3 = D[3]; N4 = D[4];}
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {NRFS = D[0];}
    void Set_External_Discretisation(std::vector<int> &D)   {NRES = D[0];}

    void Set_Dimensions(std::vector<Real> &D)               {RB = D[0];  RT = D[1]; Z = -D[2];  H1 = D[3]; R = RT;}
    void Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0];}

    //--- Geometry functions
    void Generate_Nodes();
};

}

#endif // VOLUME_OF_REVOLUTION_H
