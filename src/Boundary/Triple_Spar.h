/****************************************************************************
    BEMUse Triple spar floater geometry
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

    -> This contains the general geometry of the triple spar
    -> I will derive from tis in order to general specific geos

*****************************************************************************/


#ifndef TRIPLE_SPAR_H
#define TRIPLE_SPAR_H


//#include "Boundary_Base.h"
#include "Volume_of_Revolution.h"

namespace BEMUse
{

// Baseline Triple spar geometry

class Triple_Spar : public Boundary
{
protected:

    // Discretisation
    int NA;         // Boundary discretisation
    int NAE, NRE;   // External discretisation

    // Flags
    bool ShiftRearwards = false;

    // Geo parameters
    Real R_Leg;
    Real ScaleFac = 1.0;
    Real RFS = 1.0;
    Vector3 C1 = Vector3::Zero();   // Leg 1 origin
    Vector3 C2 = Vector3::Zero();   // Leg 2 origin
    Vector3 C3 = Vector3::Zero();   // Leg 3 origin

    // Geo Objects
    CoordSys *CSLeg1, *CSLeg2, *CSLeg3;
    Boundary *Central_Cylinder = NULL;
    Boundary *Leg1 = NULL;
    Boundary *Leg2 = NULL;
    Boundary *Leg3 = NULL;

public:

    //--- Constructor
    Triple_Spar();

//    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D);
    void Set_Auxiliary_Discretisation(std::vector<int> &D);
    void Set_External_Discretisation(std::vector<int> &D);

    void Set_Dimensions(std::vector<Real> &D);
    void Set_External_Dimensions(std::vector<Real> &D);

    void Set_Flags(std::vector<bool> &D);

    //--- Geometry functions
    void Generate_Nodes();
    void Generate_Aux_Nodes();

    void Generate_Elements();
    void Generate_Aux_Elements();

    //--- Unify multiple geos
    void Unify_Geometries();

};

class SemiSub_OC4 : public Triple_Spar
{
protected:

public:

    //--- Constructor
    SemiSub_OC4();

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D);
    void Set_Dimensions(std::vector<Real> &D);

    //--- Geometry calculations (Hard coded)
    void Set_Kin_Params();
};

}

#endif // TRIPLE_SPAR_H
