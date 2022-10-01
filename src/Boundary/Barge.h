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

    -> This is the header which contains the BEMuse Barge type

*****************************************************************************/

#ifndef BARGE_H
#define BARGE_H


#include "Boundary_Base.h"

namespace BEMUse
{

class Barge : public Boundary
{

protected:

    // Discretisation
    int NX, NY, NZ, NXFS, NYFS;

    std::vector<SP_Node> NXFace1, NXFace2, NYFace1, NYFace2, NZFace;

    // Geo vars
    Real L, W, D, LFS, WFS;

public:

    //--- Constructor
    Barge()   {}

    //--- Generate elements
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes();
    void Generate_Aux_Elements();
    void Generate_Ext_Nodes();
    void Generate_Ext_Elements();

    //--- Geometry specification
//    void Set_Discretisation(std::vector<int> &Disc)    {NX=Disc[0]; NY=Disc[1]; NZ=Disc[2]; NXFS=Disc[3]; NYFS=Disc[4];}
//    void Set_Dimensions(std::vector<Real> &Dim)       {L=Dim[0]; W=Dim[1]; D=Dim[2]; LFS=Dim[3]; WFS=Dim[4];}

    void Set_Discretisation(std::vector<int> &Disc)            {NX=Disc[0]; NY=Disc[1]; NZ=Disc[2];}
    void Set_Auxiliary_Discretisation(std::vector<int> &Disc)  {}
    void Set_External_Discretisation(std::vector<int> &Disc)   {NXFS=Disc[0]; NYFS=Disc[1];}

    void Set_Dimensions(std::vector<Real> &Dim)               {L=Dim[0]; W=Dim[1]; D=Dim[2]; }
    void Set_External_Dimensions(std::vector<Real> &Dim)      {LFS=Dim[0]; WFS=Dim[1];}

    //--- ID retrieval
    int Node_ID(int A, int Z)       {}
};


}

#endif // BARGE_H
