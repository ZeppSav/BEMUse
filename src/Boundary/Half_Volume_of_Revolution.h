/****************************************************************************
    BEMUse Half Volume of revolution
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

    -> This generates a half volume of revolution.
    -> The geometry has its axis of symmetry parallel to the free surface

*****************************************************************************/

#ifndef HALF_VOLUME_OF_REVOLUTION_H
#define HALF_VOLUME_OF_REVOLUTION_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Half_Volume_of_Revolution : public Boundary
{

protected:

    //--- Geo variables
    Real R,L;
    Real XFS = 0, YFS = 0;

    //--- Discretisation variables
    int NA, NR, NL;
    int NXAS=0, NYAS=0;
    int NXFS=0, NYFS=0;

    //--- Perimeter variables
    std::vector<Vector3> Perimeter;
    virtual void Set_Perimeter();

public:

    //--- Constructor
    Half_Volume_of_Revolution()   {}

    //--- Geometry Generation
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes();
    void Generate_Aux_Elements();
    void Generate_FreeSurface_Nodes()       {}
    void Generate_FreeSurface_Elements()    {}

    int Node_ID(int A, int Z);
    int Aux_Node_ID(int A, int Z)   {}
    int Ext_Node_ID(int A, int Z)   {}
};

}

#endif // HALF_VOLUME_OF_REVOLUTION_H
