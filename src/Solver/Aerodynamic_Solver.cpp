//---------------------------------------------------------------
//------------------Aero solver Functions------------------------
//---------------------------------------------------------------

#include "Aerodynamic_Solver.h"

namespace BEMUse
{

//--- Solver parameter specification

void Aerodynamic_Solver::Set_Parameters(std::vector<Parameter> &Params)
{
    // Parameters for the solver are specified here

    for (Parameter P : Params)
    {
        if (P.myNameis("ConstantPanels"))   PanelDist = CONSTANT;
        if (P.myNameis("BilinearPanels"))   PanelDist = BILINEAR;
    }

    // Hard-coded configurations:
    if (PanelDist==BILINEAR) std::cout << "Triangular panels must be used for bilinear strengths elements!!!" << std::endl;
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
    NPTOT = NPA;        // + NPAux; No auxiliary panels now
    NNTOT = size(Body_Nodes);
}

void Aerodynamic_Solver::Prepare_Linear_System_Const()
{
    // Prepares the linear system for the case that constant dist. panels are being used

    G = Matrix::Zero(NPTOT,NPTOT);
    H = Matrix::Zero(NPTOT,NPTOT);

    OpenMPfor
    for (int S=0; S<NPTOT; S++){
        for (int R=0; R<NPTOT; R++){
            Body_Panels[S]->Inf_SingleDoubleLayer(BC_Pos[R],G(R,S),H(R,S));
        }
    } 

    rPPLU.compute(H);
    // GMRES.compute(H);
}

void Aerodynamic_Solver::Prepare_Linear_System_Bilinear()
{
    // Prepares the linear system for the case that bilinear dist. panels are being used

    G = Matrix::Zero(NNTOT, NNTOT);
    H = Matrix::Zero(NNTOT, NNTOT);

    OpenMPfor
    for (int R = 0; R < NNTOT; R++) {
        for (int S = 0; S < NPTOT; S++) {
            Real s1, s2, s3, d1, d2, d3;
            Body_Panels[S]->Inf_SingleDoubleLayer_LinDist(BC_Pos[R], s1, s2, s3, d1, d2, d3);

            int sid1 = Body_Panels[S]->Get_Geo()->Get_Node(0)->ID;
            int sid2 = Body_Panels[S]->Get_Geo()->Get_Node(1)->ID;
            int sid3 = Body_Panels[S]->Get_Geo()->Get_Node(2)->ID;

            G(R, sid1) += s1;
            G(R, sid2) += s2;
            G(R, sid3) += s3;

            H(R, sid1) += d1;
            H(R, sid2) += d2;
            H(R, sid3) += d3;

            // Debugging flag:
            // if (std::isnan(s1) || std::isnan(s2) || std::isnan(s3) || std::isnan(d1) || std::isnan(d2) || std::isnan(d3)){
            //     Body_Panels[S]->Print_Local_Coordinates(BC_Pos[R]);
            // }
        }
    }

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Linear System diagnostics:  " << std::endl;
    std::cout << "G matrix min & max:" << G.minCoeff() csp G.maxCoeff() << std::endl;
    std::cout << "H matrix min & max:" << H.minCoeff() csp H.maxCoeff() << std::endl;

    bool GNans = G.array().isNaN().any();
    bool HNans = H.array().isNaN().any();
    if (GNans) std::cout << "G matrix contains #nans:- Check setup" << std::endl;
    if (HNans) std::cout << "H matrix contains #nans:- Check setup" << std::endl;
    if (GNans || HNans){
        std::vector<std::pair<int, int>> nan_positionsG, nan_positionsH;
        for (int i = 0; i < G.rows(); ++i) {
            for (int j = 0; j < G.cols(); ++j) {
                if (std::isnan(G(i, j))) {
                    nan_positionsG.emplace_back(i, j);
                }
                if (std::isnan(H(i, j))) {
                    nan_positionsH.emplace_back(i, j);
                }
            }
        }
    }

    if ((G.array() == 0).any()) std::cout << "G Matrix has at least one zero entry." << std::endl;
    if ((H.array() == 0).any()) std::cout << "H Matrix has at least one zero entry." << std::endl;

    std::cout << "-----------------------------" << std::endl;

    rPPLU.compute(H);    // Solution with piv LU decomposition
    // HInv = H.inverse();     // Solution with inverse inverse matrix
}

//--- Setup

void Aerodynamic_Solver::Setup(Boundary *B)
{
    // Generate surface/panels
    Create_Panels(B);

    // Specify parameters for Boundary conditions
    if (PanelDist==CONSTANT)  {for (size_t i=0; i<NPTOT; i++) BC_Nodes.push_back(Body_Panels[i]->Get_Geo()->Centroid);}
    if (PanelDist==BILINEAR)  {StdAppend(BC_Nodes,Body_Nodes);}
    Specify_BC_Positions();

    // Prepare Linear System
    if (PanelDist==CONSTANT)    Prepare_Linear_System_Const();
    if (PanelDist==BILINEAR)    Prepare_Linear_System_Bilinear();

    // Set solver parameters
//    GMRES.setMaxIterations(10);
//    GMRES.setTolerance(1.0e-10);
}

//--- Boundary conditions

void Aerodynamic_Solver::Set_External_BC(std::vector<Vector3> &Vels)
{
    // This specifies the boundary condition vector in the case that the velocity field is specified from an outside source.
    if (size(Vels)!=size(BC_Nodes)){
        std::cout << "Aerodynamic_Solver::Set_External_BC: Incorrect number of boundary condition positions." <<std::endl;
        BC.setZero();
        return;
    }

    OpenMPfor
    for (size_t i=0; i<size(BC_Nodes); i++){
        Vector3 Z = BC_Nodes[i]->Z_Axis_Global();
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
    for (int i=0; i<NPTOT; i++) Body_Panels[i]->Get_Geo()->Centroid->VWeight(0) = X(i); // Constant strength
    // BodySurface->
}

}
