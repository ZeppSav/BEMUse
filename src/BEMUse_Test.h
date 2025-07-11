/****************************************************************************
    BEMUser Testing class
    Copyright (C) 2025 Joseph Saverin j.saverin@tu-berlin.de

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

    Information on class:

    -> This class is set up to automate testing

*****************************************************************************/

#ifndef BEMUSE_TEST_H
#define BEMUSE_TEST_H

// Geometries
#include "src/Boundary/Ellipsoid.h"

// Solver types
#include "src/Solver/Surface.h"
#include "src/Solver/Aerodynamic_Solver.h"

//-----------Parent class-------------

class BEMUse_Test
{
protected:

    // Identifiers
    std::string Test_Name;
    std::vector<BEMUse::Parameter> Params;

    // Geometry
    BEMUse::Boundary *Geometry = nullptr;
    std::vector<BEMUse::SP_Node> BC_Nodes;
    std::vector<bool> SolverFlags = {true, true, true, true};   // Assume all surfaces are being generated
    BEMUse::Surface *BodySurface = nullptr;
    BEMUse::Surface *WallSurface = nullptr;
    BEMUse::Surface *SeaBedSurface = nullptr;
    BEMUse::Surface *FreeSurface = nullptr;
    std::vector<BEMUse::SP_Node> BDNodes;           int NBNodes = 0;
    std::vector<BEMUse::SP_Node> WNodes;            int NWNodes = 0;
    std::vector<BEMUse::SP_Node> SBNodes;           int NSBNodes = 0;
    std::vector<BEMUse::SP_Node> FSNodes;           int NFSNodes = 0;

    int NNodesTot = 0;

    // Unsteady tests
    int timestep = 0;
    Real t = 0.;                // Current time
    Real dt = 0.;                // Time step

    // Storage of elements
    Matrix Eta_Analytical;      // Analytical solution
    Matrix Eta_Numerical;       // BEMUse solution
    Matrix BC_vector;
    Matrix Sol_vector;

    // Solver
    // BEMUse::Rankine_Source_Solver *Solver = nullptr;
    BEMUse::Aerodynamic_Solver *Solver = nullptr;


public:

    //--- Constructor
    BEMUse_Test(std::vector<BEMUse::Parameter> &P)  {StdAppend(Params,P);}

    //--- Testing functions
    void Execute_Test()
    {
        // Set up geometry
        Geometry->Set_Parameters(Params);
        Geometry->Setup();

        // Set up solver
        Solver = new BEMUse::Aerodynamic_Solver();
        Solver->Set_Flags(SolverFlags);
        Solver->Setup(Geometry);

        // // Extract out surfaces and BC nodes
        Solver->Get_Surfaces(BodySurface, WallSurface, SeaBedSurface, FreeSurface);
        // if (BodySurface)    BodySurface->Get_Nodes(BDNodes);
        // if (WallSurface)    WallSurface->Get_Nodes(WNodes);
        // if (SeaBedSurface)  SeaBedSurface->Get_Nodes(SBNodes);
        // if (FreeSurface)    FreeSurface->Get_Nodes(FSNodes);
        // Solver->Get_BC_Nodes(BC_Nodes);
        // NNodesTot = size(BC_Nodes);

        // Execute
        Calculate();

        // Clean up
        if (Solver)     delete Solver;
        if (Geometry)   delete Geometry;
    }

    virtual void Calculate()                {}
    virtual void Output_Results()           {}
};

class Ellipsoid_Test : public BEMUse_Test
{
protected:

    Vector3 Uinf = Vector3::Zero();

public:
    //--- Constructor
    Ellipsoid_Test(std::vector<BEMUse::Parameter> &P) : BEMUse_Test(P) {
        Geometry = new BEMUse::Ellipsoid();
        SolverFlags =  {true,false,false,false};    // ONLY Ellipsoid is generated in this test case
        for (BEMUse::Parameter P : Params)
        {
            if (P.myNameis("Uinf_x"))   Uinf(0) = P.Get_Param<Real>();
            if (P.myNameis("Uinf_y"))   Uinf(1) = P.Get_Param<Real>();
            if (P.myNameis("Uinf_z"))   Uinf(2) = P.Get_Param<Real>();
        }
    }

    //--- Testing functions

    void Calculate() override
    {
        // In this test case, the boundary solution of the Ellipsoid is validated
        // for the case that the ambient inflow is along the three principle axes

        // Extract positions where boundary conditions are evaluated
        std::vector<Vector3> BCPos;
        Solver->Extract_BC_Pos(BCPos);
        NBNodes = size(BCPos);

        // Set Boundary conditions
        std::vector<Vector3> BCArray(NBNodes,Uinf);     // Generate array for constant windspeed
        Solver->Set_External_BC(BCArray);

        // Solve system
        Solver->Solve_Steady();
        std::cout << "Solved" << std::endl;
        // Solver->Get_Solution(Sol_vector);

        if (BodySurface)    BodySurface->Export_VTP();

        std::cout << "Exported" << std::endl;

        // Output solution
        // Output_Results();
    }

    void Output_Results() override
    {
        // Reference values
        int NBNodes = size(BC_Nodes);
        std::vector<Vector3> DG(NBNodes);
        Vector Xv = Sol_vector.reshaped();

        BodySurface->Calculate_PG_Coeffs();
        BodySurface->Calculate_PG_Gradients(Xv,DG);

        std::vector<Vector3> DB(NBNodes);
        BodySurface->BSpline_Coefficients();
        BodySurface->BSpline_Gradient(Xv,DB);

        // Test accuracy against analytical solution for ellipsoid
        //Analytical solution for (4,2,1) Ellipsoid
        Real alpha = 0.224700883151007738;
        Real beta = 0.569560963385673906;
        Real gamma = 1.205738153463318356;

        // Analytical solution for ellipsoid
        Vector3 gradPhi_analytical = Vector3(Uinf(0)*2./(2.-alpha), Uinf(1)*2./(2.-beta), Uinf(2)*2./(2.-gamma));

        //Error calculation
        Real L2_pressure_linear_sum = 0, L2_pressure_bspline_sum = 0, L2_velocityPotential_sum = 0;
        Real Lmax_pressure_linear = 0, Lmax_velocityPotential = 0, Lmax_pressure_bspline = 0;
        Real weight = 1.0 / NBNodes;
        Real area = Geometry->Surface_Area;
        Real L_ref = 8;  // 2 * 4
        std::vector<BEMUse::SP_Panel> Body_Panels;
        BodySurface->Get_Panels(Body_Panels);
        int NPA = size(Body_Panels);
        Real h_ref = sqrt(area / NPA);

        for (int i=0; i<NBNodes; i++) {
            // Calculate analytical solution
            Vector3 xl = BC_Nodes[i]->X_Axis_Global();
            Vector3 yl = BC_Nodes[i]->Y_Axis_Global();
            Vector3 zl = BC_Nodes[i]->Z_Axis_Global();
            Vector3 P = BC_Nodes[i]->Position_Global();

            Vector3 sol_grad_analytical = gradPhi_analytical - gradPhi_analytical.dot(zl) * zl;
            Real phi_analytical = (Uinf(0) * P[0] * alpha / (2. - alpha)) +
                                  (Uinf(1) * P[1] * beta / (2. - beta)) +
                                  (Uinf(2) * P[2] * gamma / (2. - gamma));

            // Calculate solutions due to surface interpolation (remove tangent inflow component)
            Real Un = Uinf.dot(zl);
            Vector3 Uinftang = Uinf-Un*zl;
            Real uinfx = Uinftang.dot(xl);
            Real uinfy = Uinftang.dot(yl);
            DG[i] += uinfx*xl + uinfy*yl;
            DB[i] += uinfx*xl + uinfy*yl;

            // Set visualisation
            // BC_Nodes[i]->VWeight(0) = Sol_vector(i);
            // BC_Nodes[i]->VWeight(0) = (DG[i]-sol_Analytical).norm();         // Visualise panel gradient error
            BC_Nodes[i]->VWeight(0) = (DB[i]-sol_grad_analytical).norm();    // Visualise bspline error
            // BC_Nodes[i]->VWeight(0) = (DB[i]-sol_grad_analytical).norm();    // Visualise bspline error
            // BC_Nodes[i]->VWeight(0) = (Sol_vector(i) - phi_analytical);         // Visualise panel gradient error

            Real cp_analytical = 1 - sol_grad_analytical.squaredNorm();
            Real diff_velocityPotential = abs(Sol_vector(i) - phi_analytical);
            Real diff_bspline = abs(cp_analytical - (1 - DB[i].squaredNorm()));
            Real diff_linear = abs(cp_analytical - (1 - DG[i].squaredNorm()));

            Lmax_velocityPotential = std::max(Lmax_velocityPotential, diff_velocityPotential);
            Lmax_pressure_linear = std::max(Lmax_pressure_linear, diff_linear);
            Lmax_pressure_bspline = std::max(Lmax_pressure_bspline, diff_bspline);

            L2_velocityPotential_sum += weight * pow(diff_velocityPotential, 2);
            L2_pressure_bspline_sum += weight * pow(diff_bspline, 2);
            L2_pressure_linear_sum += weight * pow(diff_linear, 2);
        }

        // Export surfaces to paraview format
        if (Uinf(0) != 0)       BodySurface->Set_Name("_Ellipsoid_Test_Inflow_X.csv");
        else if (Uinf(1) != 0)  BodySurface->Set_Name("_Ellipsoid_Test_Inflow_Y.csv");
        else if (Uinf(2) != 0)  BodySurface->Set_Name("_Ellipsoid_Test_Inflow_Z.csv");
        BodySurface->Export_VTP();
        // FreeSurface->Export_VTP();

        // Compute L2 norms
        Real L2_velocityPotential = sqrt(L2_velocityPotential_sum);
        Real errorL2_velocityPotential = L2_velocityPotential / (L_ref * Uinf.norm());
        Real errorLmax_velocityPotential = Lmax_velocityPotential / (L_ref * Uinf.norm());

        // Determine output file path
        std::string filePath;
        if (Uinf(0) != 0) filePath = "errorX.csv";
        else if (Uinf(1) != 0) filePath = "errorY.csv";
        else if (Uinf(2) != 0) filePath = "errorZ.csv";
        else {
            std::cerr << "Error: Invalid Uinf direction, file not created!" << std::endl;
            return;
        }

        std::ifstream checkFile(filePath);
        bool fileExisted = checkFile.good();
        checkFile.close();

        std::ofstream resultFile(filePath, std::ios::app);
        if (!resultFile) {
            std::cerr << "Error: Failed to open file: " << filePath << std::endl;
            return;
        }

        // Write header if file is new
        if (!fileExisted) {
            resultFile << "NPA,area,L_ref,L2_velocityPotential,errorL2_velocityPotential,L2_pressure_linear,L2_pressure_bspline,"
                          "errorLmax_velocityPotential,Lmax_pressure_bspline,Lmax_pressure_linear,h_ref\n";
        }

        resultFile << NPA << "," << area << "," << L_ref << "," << L2_velocityPotential << ","
                   << errorL2_velocityPotential << "," << std::sqrt(L2_pressure_linear_sum) << ","
                   << std::sqrt(L2_pressure_bspline_sum) << "," << errorLmax_velocityPotential << ","
                   << Lmax_pressure_bspline << "," << Lmax_pressure_linear << "," << h_ref << "\n";
    }
};

#endif // BEMUSE_TEST_H
