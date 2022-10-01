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

    -> This contains the class which imports the geometry from an stl file.

*****************************************************************************/

#ifndef WING_H
#define WING_H

#include "Boundary_Base.h"
#include "Airfoil_Profiles.h"

namespace BEMUse
{

class Wing : public Boundary
{
protected:

    // Geometry variables
    int NS, NC;
    Real Span, C;
    StdVector SpanPos, Chord, Alpha0;
    std::vector<CoordSys*> SpanwiseCS;

    //--- Geometry functions
    void SetSpan();
    void SetChord();
    void SetAlpha0();
    void SetProfiles();

    //--- Node access
    int Node_ID(int A, int Z)           {return A + (NC+1)*Z;}
    int Aux_Node_ID(int A, int Z)       {}
    int Ext_Node_ID(int A, int Z)       {}

public:

    //--- Constructor
    Wing()      {}

    //--- Geometry generation
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes()       {}
    void Generate_Aux_Elements()    {}
    void Generate_Ext_Nodes()       {}
    void Generate_Ext_Elements()    {}

    //--- Geometry specification
    void Set_Discretisation(std::vector<int> &D)            {NS = D[0]; NC = D[1];}
    void Set_Auxiliary_Discretisation(std::vector<int> &D)  {}
    void Set_External_Discretisation(std::vector<int> &D)   {}

    void Set_Dimensions(std::vector<Real> &D)               {Span = D[0]; C = D[1];}
    void Set_Auxiliary_Dimensions(std::vector<Real> &D)     {}
    void Set_External_Dimensions(std::vector<Real> &D)      {}

    void Set_Flags(std::vector<bool> &D)                    {Cosine = D[0];}
};

}

#endif // WING_H
