/****************************************************************************
    BEMUse Aerodynamic Solver
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

    -> Aerodynamic panel method solver.

*****************************************************************************/

#ifndef AERODYNAMIC_SOLVER_H
#define AERODYNAMIC_SOLVER_H

#include "Solver_Base.h"
#include "Surface.h"

namespace BEMUse
{

class Aerodynamic_Solver : public Solver
{
protected:

    //--- Problem Setup
    void Create_Panels(Boundary *B);
    void Specify_BC_Nodes();
    void Prepare_Linear_System_Const();
    void Prepare_Linear_System_Bilinear();

    //--- Problem solution
    void Set_RHS_Vec()               {}

public:

    //--- Constructor
    Aerodynamic_Solver() {}

    //--- Solver parameter specification
    void Set_Parameters(std::vector<Parameter> &Params) override;

    //--- Setup
    void Setup(Boundary *B);

    //--- Boundary conditions
    void Set_External_BC(std::vector<Vector3> &Vels)  override; // Ambient flow specified outside

    //--- Solution
    void Solve_Steady();

    //--- Post processing
    void Post_Processing(Boundary *B)       {}

    void Get_Surfaces(Surface*& S1, Surface*& S2, Surface*& S3, Surface*& S4)
    {
        S1 = BodySurface;
        // if (isWall)         S2 = WallSurface;
        // if (isSeaBed)       S3 = SeaBedSurface;
        // if (isFreeSurface)  S4 = FreeSurface;
    }

};

}

#endif // AERODYNAMIC_SOLVER_H
