//---------------------------------------------------------------
//------------------Aero solver Functions------------------------
//---------------------------------------------------------------

#include "Aerodynamic_Solver.h"

namespace BEMUse
{

//--- Problem Setup
void Aerodynamic_Solver::Create_Panels(Boundary *B)
{
    // Generate source panels.
    std::vector<SP_Geo> Geo;
    B->Get_Elements(Geo);

    NPA = Geo.size();
    for (int i=0; i<NPA; i++){
        if (Geo[i]->Get_N()==3)     Panels.push_back(std::make_shared<FlatSourceTriPanel>(Geo[i]));
        if (Geo[i]->Get_N()==4)     Panels.push_back(std::make_shared<FlatSourceQuadPanel>(Geo[i]));
    }

    std::vector<SP_Geo> Geo_Aux;
    B->Get_Aux_Elements(Geo_Aux);
    NPAux = Geo_Aux.size();
    for (int i=0; i<NPAux; i++){
        if (Geo[i]->Get_N()==3)     Panels.push_back(std::make_shared<FlatSourceTriPanel>(Geo_Aux[i]));
        if (Geo[i]->Get_N()==4)     Panels.push_back(std::make_shared<FlatSourceQuadPanel>(Geo_Aux[i]));
    }

    NPTot = NPA + NPAux;

    // Now specify panel nodes (for BC later)
    for (SP_Panel P : Panels)   Panel_Nodes.push_back(P->Get_Geo()->Centroid);

}

void Aerodynamic_Solver::Prepare_Linear_System()
{
    // Prepare the linear system

    SMat = CMatrix::Zero(NPTot,NPTot);
    DMat = CMatrix::Zero(NPTot,NPTot);

    OpenMPfor
    for (int S=0; S<NPTot; S++){
        for (int R=0; R<NPTot; R++){
            Real s, d;
            Panels[S]->Inf_SingleDoubleLayer(BC_Pos[R],s, d);
            SMat(R,S) = CReal(s,0.0);
            DMat(R,S) = CReal(d,0.0);
        }
    }
}

//--- Setup
void Aerodynamic_Solver::Setup(Boundary *B)
{
    // Prepare system
    Create_Panels(B);
    Specify_BC_Const_Pans(B);
    Prepare_Linear_System();

    // Set solver parameters
//    GMRES.setMaxIterations(10);
//    GMRES.setTolerance(1.0e-10);
}

//--- Solution
void Aerodynamic_Solver::Solve()
{
    //--- This carries out the solve for the current BC

    //--- Set BC
    RHSMat = CMatrix::Zero(NPTot,1);
    for (int i=0; i<NPTot; i++){
        Vector3 P = Panel_Nodes[i]->Position_Global();
        Vector3 N = Panel_Nodes[i]->Z_Axis_Global();
        Vector3 PCN = P.cross(N);
//        RHSMat(i) = CReal(N(2),0.0);        // Translation
        RHSMat(i) = CReal(PCN(0),0.0);      // Rotation
    }

    //--- Calc solution
    CMatrix RHSTemp = SMat*RHSMat;

    // Prepare linear system
//    PPLU.compute(DMat);                     // Prepare linear solver for solution
    GMRES.compute(DMat);

    // Solve
//    CMatrix Phi_J = PPLU.solve(RHSTemp);                       // Solution using a partial piv Lu decomposition
    CMatrix Phi_J = GMRES.solve(RHSTemp);
//    VisMat = DMat;
    std::cout << "GMRES: #iterations: " << GMRES.iterations() << ", estimated error: " << GMRES.error() << std::endl;

    RadSolArray.push_back(Phi_J);
}

}
