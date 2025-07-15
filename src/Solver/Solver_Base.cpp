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

    for (int i=0; i<NPTOT; i++){
        SP_Node C = Body_Panels[i]->Get_Geo()->Centroid;
        Vector3 P = C->Position_Global();
        Vector3 Z = C->Z_Axis_Global();
        // Real L = Panels[i]->Get_Geo()->Get_Lmax();
        Vector3 BCP = P+F*Z;

        // BC_Nodes.push_back(C);
        BC_Pos.push_back(BCP);
    }

    BC = Matrix::Zero(NPTOT,1);
}

void Solver::Specify_BC_Const_Prev(Boundary *B)
{
    // We shall assume that the pans are flat and have constant stregnths
    // This is a hacked in function just to ensure continuity between versions

    Real F = 1e-5;     // Displacement factor

    for (int i=0; i<NPTOT; i++){
        SP_Node C = Body_Panels[i]->Get_Geo()->Centroid;
        Vector3 P = C->Position_Global();
        Vector3 Z = C->Z_Axis_Global();
        // Real L = Panels[i]->Get_Geo()->Get_Lmax();
        Vector3 BCP = P+F*Z;

        BC_Nodes.push_back(C);
        BC_Pos.push_back(BCP);
    }

    // BC = Matrix::Zero(NPTOT,1);
}

void Solver::Specify_BC_Linear_Pans(Boundary *B)
{
    // We shall assume that the pans are flat and have constant stregnths

    Real F = 1e-5;     // Displacement factor

    for (int i=0; i<NNTOT; i++){
        Vector3 P = Body_Nodes[i]->Position_Global();
        Vector3 Z = Body_Nodes[i]->Z_Axis_Global();
        // Real L = Body_Panels[i]->Get_Geo()->Get_Lmax();
        Vector3 BCP = P+F*Z;
        BC_Pos.push_back(BCP);
        BC_Nodes.push_back(Body_Nodes[i]);
    }

    BC = Matrix::Zero(NNTOT,1);
}

void Solver::Specify_BC_Positions()
{
    // The positions where the boundary conditions are being evaluated are specified.
    // The factor BCF should be changed depending on where the BIE Boundary conditions are evaluted

    for (size_t i=0; i<size(BC_Nodes); i++){
        Vector3 P = BC_Nodes[i]->Position_Global();
        Vector3 Z = BC_Nodes[i]->Z_Axis_Global();
        BC_Pos.push_back(P+BCF*Z);
    }

    BC = Matrix::Zero(size(BC_Nodes),1);
}

}
