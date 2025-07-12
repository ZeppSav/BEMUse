//---------------------------------------------------------------
//------------------Aero solver Functions------------------------
//---------------------------------------------------------------

#include "Aerodynamic_Solver.h"

namespace BEMUse
{

//--- Solver parameter specification

void Aerodynamic_Solver::Set_Flags(std::vector<bool> &D)
{
    // isBody = D[0];
    // isWall = D[1];
    // isSeaBed = D[2];
    // isFreeSurface = D[3];
}

//--- Problem Setup
void Aerodynamic_Solver::Create_Panels(Boundary *B)
{
    // Generate surface representation
    std::vector<SP_Geo> BodyGeo;
    B->Get_Elements(BodyGeo);
    B->Get_Nodes(Body_Nodes);
    BodySurface = new Surface(BodyGeo,Body_Nodes,PanelDist);
    BodySurface->Set_Name("Body");
    BodySurface->Get_Panels(Body_Panels);
    NPA = size(Body_Panels);
    NPTot = NPA;        // + NPAux; No auxiliary panels now
}

void Aerodynamic_Solver::Prepare_Linear_System()
{
    // Prepare the linear system (assume for now constant panel strengths)

    G = Matrix::Zero(NPTot,NPTot);
    H = Matrix::Zero(NPTot,NPTot);

    OpenMPfor
    for (int S=0; S<NPTot; S++){
        for (int R=0; R<NPTot; R++){
            Body_Panels[S]->Inf_SingleDoubleLayer(BC_Pos[R],G(R,S),H(R,S));
        }
    }

    rPPLU.compute(H);
    // GMRES.compute(H);
}

//--- Setup
void Aerodynamic_Solver::Setup(Boundary *B)
{
    // Prepare system
    Create_Panels(B);
    if (PanelDist==CONSTANT)    Specify_BC_Const_Pans(B);
    if (PanelDist==BILINEAR)    Specify_BC_Linear_Pans(B);
    Prepare_Linear_System();

    // Set solver parameters
//    GMRES.setMaxIterations(10);
//    GMRES.setTolerance(1.0e-10);
}

//--- Boundary conditions

void Aerodynamic_Solver::Set_External_BC(std::vector<Vector3> &Vels)
{
    // This specifies the boundary condition vector in the case that the velocity field is specified from an outisde source.
    OpenMPfor
    for (int i=0; i<NPTot; i++){
        Vector3 Z = Body_Panels[i]->Get_Geo()->Centroid->Z_Axis_Global();
        BC(i) = Vels[i].dot(Z);
    }
}

//--- Solution
void Aerodynamic_Solver::Solve_Steady()
{
    //--- This carries out the solve for the previously specified BC.

    S = G*BC;               // Specify source strength
    X = rPPLU.solve(S);     // Specify perturbation potential
    // X = GMRES.solve(S);     // Specify perturbation potential

    // Store solution on panels
    for (int i=0; i<NPTot; i++) Body_Panels[i]->Get_Geo()->Centroid->VWeight(0) = X(i); // Constant strength
    // BodySurface->
}

}
