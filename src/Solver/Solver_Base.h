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
#include "../Geometry/Panel.h"
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
//#include "Eigen/unsupported/Eigen/src/IterativeSolvers/DGMRES.h"    //Note: Had to change path in Eigen folder to Eigen/Eigen/Eigenvalues

namespace BEMUse
{

class Solver
{

protected:

    //--- Solver vars
    std::vector<SP_Panel>   Panels, Panels_Aux;
    std::vector<SP_Node>    Panel_Nodes;
    StateVector             BC_Pos;
    int NPA=0, NPAux=0, NPTot=0;              // Panel counts
    int NA=0, NB=0, NR=0, NAux=0, NATot=0;

    //--- Linear system direct solver objects
    Eigen::PartialPivLU<CMatrix> PPLU;

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
    void Specify_BC_Bilin_Pans(Boundary *B);

    //--- Output files
    std::string OutputFile = "Output/BEMUse_Outputs";

public:

    //--- Constructor
    Solver() {}

    //--- Destructor
    virtual ~Solver()   {}          // Ideally nothing to clear!

    //--- Solver parameter specification
    virtual void Set_Ints(std::vector<int> &D)          {}
    virtual void Set_Real(Real D)                       {}
    virtual void Set_Reals(std::vector<Real> &D)        {}
    virtual void Set_Flags(std::vector<bool> &D)        {}
    virtual void Set_Environment(std::vector<Real> &D)  {}
    virtual void Set_OutputFilePath(std::string &D)     {OutputFile = "Output/" + D;}

    //--- Setup
    virtual void Setup(Boundary *B)                 {}

    //--- Solution
    virtual void Solve()                            {}
    virtual void Solve(Boundary *B)                 {}

    //--- Post processing
    virtual void Post_Processing(Boundary *B)       {}

    //--- Output
    virtual void Generate_Output_File(Boundary *B)  {}
    virtual void Update_Output_File()               {}

    //--- Functions for visualisation
    virtual CMatrix Get_VisMatrix()                 {}

    //--- Vars for visualisation
    std::vector<CMatrix> RadSolArray;           // Solution array
    std::vector<CMatrix> DiffSolArray;          // Solution array
    std::vector<CMatrix> AuxSolArray;           // Solution array

};

}

#endif // SOLVER_BASE_H
