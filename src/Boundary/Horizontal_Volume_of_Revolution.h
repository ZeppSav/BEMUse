/****************************************************************************
    BEMUse  Hoizontal volume of revolution
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

    -> This generates a hoizontal volume of revolution.
    -> This is useful for generating submerged cylindrical sections

*****************************************************************************/

#ifndef HORIZONTAL_VOLUME_OF_REVOLUTION_H
#define HORIZONTAL_VOLUME_OF_REVOLUTION_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Horizontal_Volume_of_Revolution : public Boundary
{

protected:

    //--- Geo variables
    Real R, Depth;
    Real L;

    //--- Perimeter variables
//    Vector Origin = Vector3(0,0,0);
    Real RFS = 0;
    int NA, NR, NL, NZ;
//    int NRFS=0, NRES=0;
    std::vector<Vector3> Perimeter;

    virtual void Create_Perimeter();

public:

    //--- Constructor
    Horizontal_Volume_of_Revolution()   {}

    //--- Geometry functions
    void Generate_Nodes();
    void Generate_Elements();
//    void Generate_Aux_Nodes();
//    void Generate_Aux_Elements();
//    void Generate_FreeSurface_Nodes();
//    void Generate_FreeSurface_Elements();

    int Node_ID(int A, int Z);
//    int Aux_Node_ID(int A, int Z)       {}
//    int Ext_Node_ID(int A, int Z);
};

class Double_Taper_Cylinder : public Horizontal_Volume_of_Revolution
{
protected:

    //--- Geo variables
    Real R1, R2;
    Real L1, L2;

    //--- Perimeter variables
    int NR, NL1, NL2;

    void Create_Perimeter();

public:

    //--- Constructor
    Double_Taper_Cylinder()   {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NA = D[0]; NR = D[1]; NL1 = D[2]; NL2 = D[3];}

    void Set_Dimensions(std::vector<Real> &D)               {R1 = D[0]; R2 = D[1]; Depth = -D[2]; L1 = D[3]; L2 = D[4];}
//    void Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0];}


};

}

#endif // HORIZONTAL_VOLUME_OF_REVOLUTION_H
