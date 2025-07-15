/****************************************************************************
    BEMUse Solver
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

    -> This is the header for the base solver class object.

*****************************************************************************/

#ifndef SOLVER_BASE_H
#define SOLVER_BASE_H

#include "../Boundary/Boundary_Base.h"
// #include "../Geometry/Panel.h"
#include "Surface.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
//#include "Eigen/unsupported/Eigen/src/IterativeSolvers/DGMRES.h"    //Note: Had to change path in Eigen folder to Eigen/Eigen/Eigenvalues

namespace BEMUse
{

class Solver
{

protected:

    //--- Solver vars
    // std::vector<SP_Panel>   Panels, Panels_Aux;
    std::vector<SP_Node>    Panel_Nodes;
    std::vector<SP_Node>    BC_Nodes;
    std::vector<Vector3>    BC_Pos;
    int NPA=0, NPAux=0;              // Panel counts
    int NA=0, NB=0, NR=0, NAux=0, NATot=0;

    DistType PanelDist = CONSTANT;          // What type of distribution is being used on the panels?

    //--- Solver vars
    Real BCF = 1.e-5;               // Boundary condition shifting factor.
    int NPTOT;                      // Total number of panels in the system
    int NNTOT;                      // Total number of surface points (on which the Green's second identitity is solved)

    //--- Solver variables ( these are updated at every timestep)
    Matrix G;                       // Matrices of influence coefficients
    Matrix H;                       // Matrices of influence coefficients
    Matrix X;                       // Solution vector
    Matrix BC;                      // Known BC for solution of linear system
    Matrix S;                       // panel strength
    Matrix HInv;                    // Inverse H matrix

    //--- Surfaces
    Surface *BodySurface;
    Surface *WallSurface;
    Surface *SeaBedSurface;
    Surface *FreeSurface;

    //--- Geometry
    std::vector<SP_Node>    Body_Nodes;     // Nodes on the solid body
    std::vector<SP_Node>    Wall_Nodes;     // Nodes of any sort of contianing wall
    std::vector<SP_Node>    SB_Nodes;       // Nodes on the sea bed (hydrodynamics)
    std::vector<SP_Node>    FS_Nodes;       // Nodes on the free surface (hydrodynamics)

    std::vector<SP_Panel>   Body_Panels;    // Panels on the solid body
    std::vector<SP_Panel>   Wall_Panels;    // Panels of any sort of contianing wall
    std::vector<SP_Panel>   SB_Panels;      // Panels on the sea bed (hydrodynamics)
    std::vector<SP_Panel>   FS_Panels;      // Panels on the free surface (hydrodynamics)

    //--- Linear system direct solver objects
    Eigen::PartialPivLU<CMatrix> PPLU;
    Eigen::PartialPivLU<Matrix> rPPLU;

    //--- Linear system iterative solver objects
//    Eigen::ConjugateGradient<CMatrix, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> CG;
//    Eigen::BiCGSTAB<CMatrix, Eigen::IdentityPreconditioner> BICG;
//    Eigen::DGMRES<CMatrix, Eigen::DiagonalPreconditioner<CMatrix>> DGMRES;             // Iterative solver (gmres)  << gmres.iterations() << gmres.error()
//    Eigen::MINRES<CMatrix> MINRES;
//    Eigen::DGMRES<CMatrix> DGMRES;
    Eigen::GMRES<CMatrix> GMRES;

    //--- Problem Setup
    virtual void Create_Panels(Boundary *B)         {}
//    virtual void Specify_BC_Positions(Boundary *B)  {}
    virtual void Prepare_Linear_System() {}

    //--- Create boundary condition points
    void Specify_BC_Const_Pans(Boundary *B);
    void Specify_BC_Const_Prev(Boundary *B);
    void Specify_BC_Linear_Pans(Boundary *B);
    virtual void Specify_BC_Positions();

    //--- Output files
    std::string OutputDirectory = "Output";
    std::string OutputPath = "Output/BEMUse_Outputs";

public:

    //--- Constructor
    Solver() {}

    //--- Destructor
    virtual ~Solver()   {}          // Ideally nothing to clear!

    //--- Solver parameter specification
    virtual void Set_Parameters(std::vector<Parameter> &Params) {}
    virtual void Set_Real(Real D)                               {}
    virtual void Set_OutputFilePath(std::string &D)     {OutputPath = OutputDirectory + "/" + D;}

    //--- Setup
    virtual void Setup(Boundary *B)                 {}

    //--- Boundary conditions
    virtual void Set_BC(Matrix &BCin)           {BC  = BCin;}
    virtual void Set_External_BC(std::vector<Vector3> &Vels)    {}  // Ambient flow specified outside
    virtual void Extract_BC_Pos(std::vector<Vector3> &Pos)      {StdAppend(Pos,BC_Pos);}

    //--- Solution
    virtual void Solve()                            {}
    virtual void Solve_Steady()                     {}
    virtual void Solve(Boundary *B)                 {}

    //--- Post processing
    virtual void Post_Processing(Boundary *B)       {}

    //--- Output
    virtual void Generate_Output_File(Boundary *B)  {}
    virtual void Update_Output_File()               {}

    //--- Functions for visualisation
    virtual CMatrix Get_VisMatrix()                 {}

    //--- Vars for visualisation
    std::vector<CMatrix> RadSolArray;           // Solution array radiation potential
    std::vector<CMatrix> DiffSolArray;          // Solution array scattering potential
    std::vector<CMatrix> FS_Rad_SolArray;       // Solution array Free surface elevation due to radiation
    std::vector<CMatrix> FS_Scat_SolArray;      // Solution array Free surface elevation due to scattering

};

}

#endif // SOLVER_BASE_H
