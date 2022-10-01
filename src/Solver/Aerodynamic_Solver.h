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

namespace BEMUse
{

class Aerodynamic_Solver : public Solver
{
protected:

    //--- Linear system arrays
    CMatrix SMat;    // Influence coefficient matrix (source terms)
    CMatrix DMat;    // Influence coefficient matrix (source terms)
    CMatrix RHSMat;  // RHS vector which contains the Boundary condition to be solved for
    CMatrix SolMat;  // Solution vector

    //--- Problem Setup
    void Create_Panels(Boundary *B);
    void Specify_BC_Nodes();
    void Prepare_Linear_System();

    //--- Problem solution
    void Set_RHS_Vec()               {}

public:

    //--- Constructor
    Aerodynamic_Solver() {}

    //--- Setup
    void Setup(Boundary *B);

    //--- Solution
    void Solve();

    //--- Post processing
    void Post_Processing(Boundary *B)       {}

};

}

#endif // AERODYNAMIC_SOLVER_H
