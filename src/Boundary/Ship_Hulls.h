/****************************************************************************
    BEMUse Ship hulls
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

    -> This generates a volume of revolution. The variable "Open" specifies whether the
    -> Geometry is a closed volume of revolution or not.

*****************************************************************************/

#ifndef SHIP_HULLS_H
#define SHIP_HULLS_H

#include "Boundary_Base.h"

namespace BEMUse
{

class Ship_Hull : public Boundary
{
protected:

    // Geo Params
    int NX, NZ;
    Real L=10., B=1., T=1.;

    // Free surface vars
    Real XFSMax=0, XFSMin=0;    // X FS grid max, X FS grid min
    Real XHS=0, XHK=0;          // X hull starboard, X hull keel
    Real YFS=0;                 // The Y value of the free surface grid (assume constant for now)
    std::vector<Vector3> Hull_Perim, TriHull_Perim;    // Positions of starting x & y-points for free surface mesh
    int NXK=0, NXH=0, NXS=0, NYFS=0;

    std::vector<SP_Node> LeftNodes, RightNodes;
    std::vector<SP_Node> TriLeftNodes, TriRightNodes;
    std::vector<SP_Node> Tri_FS_Nodes;

public:

    //--- Constructor
    Ship_Hull()       {}

    //--- Geometry functions
    virtual void Generate_Nodes()                   {}
    virtual void Generate_Aux_Nodes()               {}
    virtual void Generate_Elements();
    virtual void Generate_FreeSurface_Elements()    {}

    //--- Individual functions for quadratic and triangular elements
    void Generate_Quadratic_Elements();
    void Generate_Triangular_Elements();
};

// Possible alternatives:
// S175, Kriso Container, Hamburg Test Case (HTC)

class Wigley_Hull : public Ship_Hull
{

public:

    //--- Constructor
    Wigley_Hull()       {}

    //--- Geometry functions
    void Set_Parameters(std::vector<Parameter> &Params) override;

    void Generate_Nodes();
    void Generate_Aux_Nodes();
    void Generate_FreeSurface_Nodes();

    void Generate_Aux_Elements();
    void Generate_FreeSurface_Elements();
    void Correct_Panel_Orientation(std::vector<SP_Geo> *Panels, std::vector<SP_Node> *Nds);

    void Hull_Node(Real xn, Real yn, Vector3 &p, Vector3 &n);
    void FS_Node(Real xn, Real yn, Vector3 &p, int zone);

};

}

#endif // SHIP_HULLS_H
