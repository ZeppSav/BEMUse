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
    void Generate_FreeSurface_Nodes();
    void Generate_FreeSurface_Elements();

    //--- Geometry specification
    void Set_Parameters(std::vector<Parameter> &Params) override;

    //--- ID retrieval
    int Node_ID(int A, int Z)       {}
};


}

#endif // BARGE_H
