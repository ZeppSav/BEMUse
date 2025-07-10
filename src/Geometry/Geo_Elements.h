/****************************************************************************
    Node Class
    Copyright (C) 2021 Joseph Saverin j.saverin@tu-berlin.de

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

    -> This is the implementation of the node type, which inherits from the coordinate system

*****************************************************************************/

#ifndef GEO_ELEMENTS_H
#define GEO_ELEMENTS_H

#include "Node.h"

namespace BEMUse
{

class Geometry_Element
{

protected:


public:

    //--- Constructors
    Geometry_Element()      {}

    //--- Geometry construction
    virtual void Set_Centroid()         {}
    virtual void Set_Quad_Nodes()       {}
    virtual void Reorder_Nodes()        {}

    //--- Geometry construction
    std::vector<SP_Node>    Nodes;
    std::vector<SP_Node>    QNodes;
    SP_Node Centroid = nullptr;

    //--- Getters
    int     Get_N()     {return Nodes.size();}
    Real    Get_Lmin();
    Real    Get_Lmax();
    SP_Node Get_Node(int i) {return Nodes[i];}

    //--- Public data
    Real Area = 0.0;

//    //--- Specification data
//    bool isDipole = false;

    //--- Visualisation
    void    Get_Nodal_Positions(StateVector &S) {for (SP_Node N : Nodes) S.push_back(N->Position_Global());}

    StateVector Get_Nodal_Positions()
    {
        StateVector S;
        for (SP_Node N : Nodes) S.push_back(N->Position_Global());
        return S;
    }
};

class Line_Element : public Geometry_Element
{

    Real dL;
    SP_Node NodeA;
    SP_Node NodeB;

public:

    //--- Constructors
    Line_Element(SP_Node Node1, SP_Node Node2)  {Nodes.push_back(Node1); Nodes.push_back(Node2);}

};

class Area_Element : public Geometry_Element
{
protected:

    Real dA;

public:

    //--- Constructors
    Area_Element()  {}

    Vector3 Get_normal()    {return Vector3::Zero();}
};

class Tri_Element : public Area_Element
{
protected:

    int N=2;

public:

    //--- Constructors
    Tri_Element(SP_Node N1, SP_Node N2, SP_Node N3)
    {
        Nodes.push_back(N1);
        Nodes.push_back(N2);
        Nodes.push_back(N3);
    }

    //--- Geometry construction
    void Set_Centroid();
    void Set_Quad_Nodes();
    void Reorder_Nodes() {
        std::swap(Nodes[0],Nodes[2]);
    }
};

class Quad_Element : public Area_Element
{
protected:

    //--- Geometry construction
    Cart_ID2 Set_Quad_Shape();
    int N=2;

public:

    //--- Constructors
    Quad_Element(SP_Node N1, SP_Node N2, SP_Node N3, SP_Node N4)
    {
        Nodes.push_back(N1);
        Nodes.push_back(N2);
        Nodes.push_back(N3);
        Nodes.push_back(N4);
    }

    //--- Geometry construction
    void Set_Centroid();
    void Set_Quad_Nodes();
    void Reorder_Nodes() {
        std::swap(Nodes[0],Nodes[3]);
        std::swap(Nodes[1],Nodes[2]);
    }
};

typedef std::shared_ptr<Geometry_Element> SP_Geo;

}

#endif // GEO_ELEMENTS_H
