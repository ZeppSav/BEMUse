/****************************************************************************
    BEMUse Solver
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

    Information on file:

    -> The surface class represents discrete surfaces in BEMUse and contains the funcionality to
    -> e.g. find gradients over a surface.

*****************************************************************************/


#ifndef SURFACE_H
#define SURFACE_H

#include "Solver_Base.h"

namespace BEMUse
{

class Surface
{
    //--- Geometric Objects
    int NP = 0;                     // Number of panels
    int NN = 0;                     // Number of Nodes
    std::vector<SP_Node> Nodes;
    std::vector<SP_Panel> Panels;

    // Surface distribution: Panel gradient approach
    std::vector<std::vector<int>> PG_Conn;          // Connectivity IDs
    std::vector<std::vector<Real>> PG_Angles;       // Connectivity angles
    std::vector<Matrix> PG_Mats;                    // Matrices for the surface interpolation

    // Surface distribution: B-Spline approach
    std::vector<std::vector<int>> BSP_Conn;         // Neighbor IDs
    std::vector<std::vector<int>> BSP_Conn2;        // Second neighbor IDs
    std::vector<Matrix> BSP_Mats;                   // Matrices for the interpolation
    std::vector<Matrix> BSP_RBF_Mats;               // Matrices for the specification of the gradient vectors

    // Visualisation
    std::string OutputDirectory = "Output";
    std::string OutputPath = "Output/BEMUse_Outputs";
    std::string SurfaceName;

public:
    // Constructor
    Surface(std::vector<SP_Geo> &SurfPans, std::vector<SP_Node> &SurfNodes);

    // Surface distribution: Panel gradient approach
    void Specify_Connectivity_Panel_Gradient();
    void Calculate_PG_Coeffs();
    void Calculate_PG_Gradients(Vector &Field, std::vector<Vector3> &Gradient);

    // Surface distribution: B-Spline approach
    void Specify_Connectivity_BSpline();
    void BSpline_Coefficients();
    void BSpline_Gradient(Vector &Field, std::vector<Vector3> &Gradients);

    // Surface integration
    Real Integrate_Quantities_Tri(Matrix &Field, Real exp);

    // Timestepping vars
    int currentTimeStep=0;

    // Flow quantities
    Matrix purtEta;

    // Setters & Getters
    void Get_Panels(std::vector<SP_Panel> &Pans)    {StdAppend(Pans,Panels);}
    void Get_Nodes(std::vector<SP_Node> &Nds)       {StdAppend(Nds,Nodes);}
    void Set_Name(std::string Name)                 {SurfaceName = Name;}

    // Visualisation
    bool OutputWaveHeight = false;
    void Export_VTP();
};


}

#endif // SURFACE_H
