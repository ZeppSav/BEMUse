/****************************************************************************
    BEMUse Ship hulls
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

#ifndef SHIP_HULLS_H
#define SHIP_HULLS_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Ship_Hull : public Boundary
{
protected:

    // Geo Params
    int NX, NZ;
    Real L, B, T;


public:

    //--- Constructor
    Ship_Hull()       {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NX = D[0]; NZ = D[1];}
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {}
    void Set_External_Discretisation(std::vector<int> &D)   {}

    void Set_Dimensions(std::vector<Real> &D)               {L = D[0]; B = D[1]; T = D[2];}
    void Set_External_Dimensions(std::vector<Real> &D)      {}

    //--- Geometry functions
    virtual void Generate_Nodes()               {}
    virtual void Generate_Aux_Nodes()           {}
    virtual void Generate_Elements();
    virtual void Generate_Aux_Elements()        {}

    //--- ID retrieval
    int Node_ID(int X, int Z, int L);
};

// Possible alternatives:
// S175, Kriso Container,  Hamburg Test Case (HTC)

class Wigley_Hull : public Ship_Hull
{


public:

    //--- Constructor
    Wigley_Hull()       {}

    //--- Geometry functions
    void Generate_Nodes();

};

}

#endif // SHIP_HULLS_H
