/****************************************************************************
    BEMUse Ellipsoid
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

    -> This generates an ellipsoid volume

*****************************************************************************/

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Ellipsoid : public Boundary
{
protected:

    // Discretisation
    int NA, NZ;

    // Geo vars
    Real a,b,c,Depth;

    // Analytical solution
    Real alpha0, beta0, gamma0;

public:

    //--- Constructor
    Ellipsoid()   {}

    //--- Geometry generation
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes()       {}
    void Generate_Aux_Elements()    {}
    void Generate_Ext_Nodes()       {}
    void Generate_Ext_Elements()    {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NA = D[0]; NZ = D[1];}
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {}
    void Set_External_Discretisation(std::vector<int> &D)   {}

    void Set_Dimensions(std::vector<Real> &D)               {a = D[0]; b = D[1]; c = D[2]; Depth = -D[3];}
    void Set_External_Dimensions(std::vector<Real> &D)      {}

    //--- ID retrieval
    int Node_ID(int A, int Z);
//    int Aux_Node_ID(int A, int Z)       {}
};

class Semi_Ellipsoid : public Boundary      // Eg. Hemisphere
{

protected:

    // Discretisation
    int NA, NZ;
    int NRFS = 0;
    int NRES=0, NAES=0;
    // Geo vars
    Real a,b,c, RFS;

public:

    //--- Constructor
    Semi_Ellipsoid()   {}

    //--- Generate elements
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes();
    void Generate_Aux_Elements();
    void Generate_Ext_Nodes();
    void Generate_Ext_Elements();

    //--- Geometry specification
//    void Set_Discretisation(std::vector<int> &D)    {NA = D[0]; NZ = D[1]; NRFS = D[2];}
//    void Set_Dimensions(std::vector<Real> &D)       {a = D[0]; b = D[1]; c = D[2];}

    void Set_Discretisation(std::vector<int> &D)            {NA = D[0]; NZ = D[1];}
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {NRFS = D[0];}
    void Set_External_Discretisation(std::vector<int> &D)   {NAES = D[0]; NRES = D[1]; }

    void Set_Dimensions(std::vector<Real> &D)               {a = D[0]; b = D[1]; c = D[2];}
    void Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0];}

    //--- ID retrieval
    int Node_ID(int A, int Z);
    int Ext_Node_ID(int A, int Z);
//    int Aux_Node_ID(int A, int Z)       {}
};

}

#endif // ELLIPSOID_H
