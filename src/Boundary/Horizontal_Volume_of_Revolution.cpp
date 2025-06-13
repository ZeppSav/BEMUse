//---------------------------------------------------------------
//------------Horizontal VOR Functions---------------------------
//---------------------------------------------------------------

#include "Horizontal_Volume_of_Revolution.h"

namespace BEMUse
{

//---------------------------------
// Horizontal_Volume_of_Revolution
//---------------------------------

//--- Geometry functions

void Horizontal_Volume_of_Revolution::Generate_Nodes()
{
    // This function generates the nodes which shall be used to generate
    // the geometry elements for an open volume of revolution.
    // The axis of symmetry is the z-axis

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    //---- Generate perimeter to be revolved
    Create_Perimeter();

    //---- Generate coordinates

    NZ = Perimeter.size();

    Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(-0.5*L,0,Depth)));   // End node 1

    for (int i=1; i<NZ-1; i++)          //
    {
        int NAR = NA;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;
            Real x = Perimeter[i](0), r = Perimeter[i](2);
            Vector3 BP  ( x, r*cos(u), r*sin(u) + Depth);
            Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(0.5*L,0,Depth)));     // End node 2

    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;
}

void Horizontal_Volume_of_Revolution::Create_Perimeter()
{
    // Template class simply creates a cylinder
    if (Cosine){
        for (int i=0; i<NR+1; i++)  {Real f = PiCosFac(i,NR+1);     Perimeter.push_back(Vector3(-0.5*L,0,f*R));}
        for (int i=1; i<NL; i++)    {Real f = PiCosFac(i,NL+1);       Perimeter.push_back(Vector3(-0.5*L+f*L,0,R));}
        for (int i=0; i<NR+1; i++)  {Real f = PiCosFac(i,NR+1);     Perimeter.push_back(Vector3(0.5*L,0,R*(1.0-f)));}
    }
    else{
        for (int i=0; i<NR+1; i++)  {Real f = 1.0*i/NR;             Perimeter.push_back(Vector3(-0.5*L,0,f*R));}
        for (int i=1; i<NL; i++)    {Real f = 1.0*i/NL;             Perimeter.push_back(Vector3(-0.5*L+f*L,0,R));}
        for (int i=0; i<NR+1; i++)  {Real f = 1.0*i/NR;             Perimeter.push_back(Vector3(0.5*L,0.0,R*(1.0-f)));}
    }

}

void Horizontal_Volume_of_Revolution::Generate_Elements()
{
    // This function generates the elements for an open volume of revolution

    // End 1 (triangular elements)
    for (int i=0; i<NA; i++)        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                                                        Nodes[Node_ID(i+1,1)],
                                                                                        Nodes[Node_ID(i,1)]));

    // Central quadratic elements
    for (int z=1; z<NZ-2; z++){
        for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,z)],
                                                                                        Nodes[Node_ID(i+1,z)],
                                                                                        Nodes[Node_ID(i+1,z+1)],
                                                                                        Nodes[Node_ID(i,z+1)]));
    }

    // End 2 (triangular elements)
    for (int i=0; i<NA; i++)        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,NZ)],
                                                                                        Nodes[Node_ID(i,NZ-2)],
                                                                                        Nodes[Node_ID(i+1,NZ-2)]));
}

int Horizontal_Volume_of_Revolution::Node_ID(int A, int Z)
{
    // Returns the node ID based on the azimuthal and axial position
    if (Z==0)           return 0;
    else if (Z==NZ)     return NA*(NZ-2)+1;
    else {
        if (A>=NA)      return 1+(Z-1)*NA+(A-NA);
        else if (A<0)   return 1+(Z-1)*NA+(A+NA);
        else            return 1+(Z-1)*NA+A;
    }
}

//int Horizontal_Volume_of_Revolution::Ext_Node_ID(int A, int Z)
//{
//    // Returns the node ID based on the azimuthal and axial position
//    if (A>=NA)      return Z*NA + (A-NA);
//    else            return Z*NA + A;
//}

//---------------------------------
// Double Taper Cylinder
//---------------------------------

void Double_Taper_Cylinder::Create_Perimeter()
{
    // This creates the tapering cylinder geometry.
    // Note:    L2 is the distance between taper points
    //          L1 is the distance between the two ends
    Vector3 P1(-0.5*L1,0,0);
    Vector3 P2(-0.5*L1,0,R1);
    Vector3 P3(-0.5*L2,0,R2);
    Vector3 P4( 0.5*L2,0,R2);
    Vector3 P5( 0.5*L1,0,R1);
    Vector3 P6( 0.5*L1,0,0);

    if (Cosine){
//        for (int i=0; i<NR+1; i++)  {Real f = PiCosFac(i,NR+1);     Perimeter.push_back(Vector3(-0.5*L,0,f*R));}
//        for (int i=1; i<NL; i++)    {Real f = PiCosFac(i,NL+1);       Perimeter.push_back(Vector3(-0.5*L+f*L,0,R));}
//        for (int i=0; i<NR+1; i++)  {Real f = PiCosFac(i,NR+1);     Perimeter.push_back(Vector3(0.5*L,0,R*(1.0-f)));}
        for (int i=0; i<NR+1; i++)  {Real f = PiCosFac(i,NR+1);     Perimeter.push_back(P1 + f*(P2-P1));}
        for (int i=1; i<NL1+1; i++) {Real f = PiCosFac(i,NL1+1);     Perimeter.push_back(P2 + f*(P3-P2));}
        for (int i=1; i<NL2+1; i++) {Real f = PiCosFac(i,NL2+1);     Perimeter.push_back(P3 + f*(P4-P3));}
        for (int i=1; i<NL1+1; i++) {Real f = PiCosFac(i,NL1+1);     Perimeter.push_back(P4 + f*(P5-P4));}
        for (int i=1; i<NR+1; i++)  {Real f = PiCosFac(i,NR+1);     Perimeter.push_back(P5 + f*(P6-P5));}
    }
    else{
        for (int i=0; i<NR+1; i++)  {Real f = 1.0*i/NR;         Perimeter.push_back(P1 + f*(P2-P1));}
        for (int i=1; i<NL1+1; i++) {Real f = 1.0*i/NL1;        Perimeter.push_back(P2 + f*(P3-P2));}
        for (int i=1; i<NL2+1; i++) {Real f = 1.0*i/NL2;        Perimeter.push_back(P3 + f*(P4-P3));}
        for (int i=1; i<NL1+1; i++) {Real f = 1.0*i/NL1;        Perimeter.push_back(P4 + f*(P5-P4));}
        for (int i=1; i<NR+1; i++)  {Real f = 1.0*i/NR;         Perimeter.push_back(P5 + f*(P6-P5));}
    }

    L = L1;

}

}
