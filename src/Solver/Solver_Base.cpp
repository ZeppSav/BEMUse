//---------------------------------------------------------------
//------------------Solver Base Functions------------------------
//---------------------------------------------------------------

#include "Solver_Base.h"

namespace BEMUse
{

//--- Create boundary condition points

void Solver::Specify_BC_Const_Pans(Boundary *B)
{
    // We shall assume that the pans are flat and have constant stregnths

    Real F = 1e-5;     // Displacement factor

    for (int i=0; i<NPTot; i++){
        Vector3 P = Panels[i]->Get_Geo()->Centroid->Position_Global();
        Vector3 Z = Panels[i]->Get_Geo()->Centroid->Z_Axis_Global();
        Real L = Panels[i]->Get_Geo()->Get_Lmax();
        Vector3 BCP = P+F*Z;
        BC_Pos.push_back(BCP);
    }
}

}
