/****************************************************************************
    BEMUse Hydrodynamic Solver
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

    -> Hydodynamic panel method solver.

*****************************************************************************/

#ifndef HYDRODYNAMIC_SOLVER_H
#define HYDRODYNAMIC_SOLVER_H

#include "Solver_Base.h"

namespace BEMUse
{

enum WATER_DEPTH {Infinite, Finite};

class Hydrodynamic_Radiation_Solver : public Solver
{
    //--- Linear system arrays
    CMatrix SSrcMat, SSrcReflMat, SWaveMat;     // Influence coefficient matrix (source terms)
    CMatrix DSrcMat, DSrcReflMat, DWaveMat;     // Influence coefficient matrix (dipole terms)
    CMatrix ReflNormMat;                        // Normal of the reflected matrix
    CMatrix N_k_Mat, DPhi_J_DN;     // RHS vector which contains the Boundary conditions for the radiation problem
    CMatrix VisMat;
    Matrix PanArea;                             // Area of panels

    //--- Kochin functions
    CMatrix KochinRad, KochinDiff;

    //--- Free surface elevation
    std::vector<SP_Node> Wave_Nodes;
    CMatrix SMatExt, SMatReflExt, SWaveMatExt;
    CMatrix DMatExt, DMatReflExt, DWaveMatExt;
    CMatrix ExtRadMat, ExtDiffMat;  // Solution vector

    //--- Solution matrices
    Matrix A_ij = Matrix(6,6);
    Matrix B_ij = Matrix(6,6);
    Matrix B_ijOm = Matrix(6,6);
    CMatrix FK_i;                   // Froude-Krylov forces
    CMatrix SC_i;                   // Scattering forces
    CMatrix F_i;                    // Excitation forces
    CMatrix D_i;                    // Diffraction forces
    CMatrix RAO_i;                  // Array of RAO s

    //--- Post processing matrices
    CMatrix AddedMassMat;                       // Matrix to calculate the added mass / damping terms
    CMatrix Phi_I, DPhi_I_DN;                   // Incident potential + gradient terms
    CMatrix KochinE, KochindEdn;                // Kochin terms azimtuhal integration
    CMatrix KochinPiE, KochinPidEdn;            // Kochin terms const pi

    //--- Panels
    std::vector<SP_Panel> Source_Panels, Refl_Source_Panels, Wave_Panels;
    std::vector<SP_Geo> FS_Geo;

    //--- Problem Setup
    void Create_Panels(Boundary *B);
    void Specify_BC_Nodes();
    void Prepare_Linear_System();
    void Prepare_Linear_System_Wave_Terms();
    void Prepare_FS_Linear_System(Boundary *B);
    void Prepare_FS_Linear_System_Wave_Terms();
    void Prepare_PostProcessing_Mats();

    void Calculate_Mass_Matrix(Boundary *B);
    void Calculate_Hydrostatic_Stiffness_Matrix(Boundary *B);

    //--- Problem solution
    void Set_RHS_Mat();

    //--- Processing
    void Set_Incident_Potential_Mats();
    void Calc_Added_Mass(int DOF_Analysis);
    void Calc_Excitation_Force(int NB);
    void Set_Kochin_Mats();
    void Calc_Phi_Incident(Real &Beta, Vector3 &P_Glob, CReal &Phi_i, CReal &dPhi_iX,  CReal &dPhi_iY,  CReal &dPhi_iZ);
    Real Calc_Root_Finite_Depth(Real &O, Real &H);

    //--- Parameters
    Real Kappa = 0.0;           // Wavenumber
    Real Frequency = 0.0;       // Frequency of oscillation
    Real Period = 0.0;
    Real Omega = 0.0;
    WATER_DEPTH Depth = Infinite;
    Real H = 10000;             // Water depth
    int t_DOF, NDOF = 6;        // Which degree of free is oscillating now?
    StdVector Frequency_List;
    bool IFR = false;
    Real Rho = Rho_wat;         // Density of fluid
    Real Grav = Gravity;        // Density of fluid

    //--- Kinematic Parameters
    Matrix M_ij = Matrix::Zero(6,6);
    Matrix C_ij = Matrix::Zero(6,6), C_ij_NM = Matrix::Zero(6,6);

    //--- Incident wave parameters
    Real Beta = 0;          // What is the current incoming wave angle?
    int NBeta = 0;          // How many incident wave angles are we investigating?
    StdVector BetaArray;    // Array of incoming wave angles
    int NKoch = 36;         // Which angular discretisation are we using to calculate the Kochin functions?

    //--- Output parameters
    Real SolEps = 1e-5;
    Real FreqInf = -1.0;

    //--- Output functions
    void Generate_Output_File_BEMUse(Boundary *B);
    void Generate_Output_File_WAMIT(Boundary *B);
    void Export_Wave_Height();
    void Update_Output_File_BEMUse();
    void Update_Output_File_WAMIT();

public:

    //--- Constructor
    Hydrodynamic_Radiation_Solver() {}

    //--- Setup
    void Setup(Boundary *B);

    //--- Solver parameter specification
    void Set_Flags(std::vector<bool> &D)        {IFR = D[0];}
    void Set_Real(Real  D)                      {Frequency = D;}
    void Set_Reals(std::vector<Real> &D)        {StdAppend(BetaArray,D);}
    void Set_Ints(std::vector<int> &D)          {NKoch = D[0];}
    void Set_Environment(std::vector<Real> &D)  {Rho = D[0]; Grav = D[1];}

    //--- Solution
    void Solve();

    //--- Post processing
    void Post_Processing(Boundary *B)       {}

    //--- Output
    void Generate_Output_File(Boundary *B);
    void Update_Output_File();

    //--- Functions for visualisation
    CMatrix Get_VisMatrix()                        {return VisMat;}
};

}

#endif // HYDRODYNAMIC_SOLVER_H
