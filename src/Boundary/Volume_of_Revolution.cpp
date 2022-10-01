//---------------------------------------------------------------
//------------------ VOR Functions-------------------------------
//---------------------------------------------------------------

#include "Volume_of_Revolution.h"

namespace BEMUse
{

//--- Geometry functions

void Volume_of_Revolution::Generate_Nodes()
{
    // This function generates the nodes which shall be used to generate
    // the geometry elements for an open volume of revolution.
    // The axis of symmetry is the z-axis

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    //---- Generate coordinates

    NZ = Perimeter.size()-1;

    for (int i=0; i<NZ+1; i++)          //
    {
        int NAR = NA;
        if (i==0)   NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;
            Real r = Perimeter[i](0), z = Perimeter[i](2);
//            Vector3 BP  ( r*cos(u) + Origin(0), r*sin(u) + Origin(1), z + Origin(2));
            Vector3 BP  ( r*cos(u), r*sin(u), z);
            Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

}

void Volume_of_Revolution::Generate_Aux_Nodes()
{
    // This function generates the nodes necessary for the interior free surface panels
    if (NRFS==0) return;

    for (int i=0; i<NRFS+1; i++)
    {
        Real r;
        if (Cosine) r = R*PiCosFac(i,NRFS+1);   // Cosine
        else        r = R*i/NRFS;           // Linear

        int NAR = NA;
        if (i==0)   NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;
//            Vector3 BP  ( r*cos(u) + Origin(0), r*sin(u) + Origin(1), Origin(2));
            Vector3 BP  ( r*cos(u), r*sin(u),0);
            Aux_Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    for (int i=0; i<Aux_Nodes.size(); i++)  Aux_Nodes[i]->ID = i;

//    std::cout << "Aux_Nodes Size: " <<  Aux_Nodes.size() << std::endl;
}

void Volume_of_Revolution::Generate_Ext_Nodes()
{
    // A volume of revolution has as it's interesting plane a circle.

//    NRES = NRFS;  // Hack!

    for (int i=0; i<NRES+1; i++)
    {
        Real r = R*1.001+RFS*R*i/NRES;      // Linear
//        Real r = RFS*PiCosFac(i,NRES+1);   // Cosine
        for (int j = 0; j<NA; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NA;
//            Vector3 BP  ( r*cos(u) + Origin(0), r*sin(u) + Origin(1), Origin(2));
            Vector3 BP  ( r*cos(u), r*sin(u), 0);
            Ext_Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    for (int i=0; i<Ext_Nodes.size(); i++) Ext_Nodes[i]->ID = i;
}

void Volume_of_Revolution::Generate_Elements()
{
    // This function generates the elements for an open volume of revolution

    // Base (triangular elements)
    for (int i=0; i<NA; i++)        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                                                        Nodes[Node_ID(i+1,1)],
                                                                                        Nodes[Node_ID(i,1)]));

    // Central quadratic elements
    for (int z=1; z<NZ; z++){
        for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,z)],
                                                                                        Nodes[Node_ID(i+1,z)],
                                                                                        Nodes[Node_ID(i+1,z+1)],
                                                                                        Nodes[Node_ID(i,z+1)]));
}

}

void Volume_of_Revolution::Generate_Aux_Elements()
{
    // This function generates the elements for an open volume of revolution
    if (NRFS==0) return;

    // Base (triangular elements)
    for (int i=0; i<NA; i++)        Aux_Elements.push_back(std::make_shared<Tri_Element>(   Aux_Nodes[Node_ID(0,0)],
                                                                                            Aux_Nodes[Node_ID(i,1)],
                                                                                            Aux_Nodes[Node_ID(i+1,1)]));

    // Central quadratic elements
    for (int z=1; z<NRFS; z++){
        for (int i=0; i<NA; i++)    Aux_Elements.push_back(std::make_shared<Quad_Element>(  Aux_Nodes[Node_ID(i,z+1)],
                                                                                            Aux_Nodes[Node_ID(i+1,z+1)],
                                                                                            Aux_Nodes[Node_ID(i+1,z)],
                                                                                            Aux_Nodes[Node_ID(i,z)]));
    }

    for (int i=0; i<Aux_Elements.size(); i++)  Aux_Elements[i]->Set_Centroid();
}

void Volume_of_Revolution::Generate_Ext_Elements()
{
    // Central quadratic elements
    for (int z=0; z<NRES; z++){
        for (int i=0; i<NA; i++)    Ext_Elements.push_back(std::make_shared<Quad_Element>(  Ext_Nodes[Ext_Node_ID(i,z+1)],
                                                                                            Ext_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                                            Ext_Nodes[Ext_Node_ID(i+1,z)],
                                                                                            Ext_Nodes[Ext_Node_ID(i,z)]));
    }

    for (int i=0; i<Ext_Elements.size(); i++)  Ext_Elements[i]->Set_Centroid();
}

int Volume_of_Revolution::Node_ID(int A, int Z)
{
    // Returns the node ID based on the azimuthal and axial position
    if (Z==0) return 0;
    else {
        if (A>=NA)      return 1+(Z-1)*NA+(A-NA);
        else if (A<0)   return 1+(Z-1)*NA+(A+NA);
        else            return 1+(Z-1)*NA+A;
    }
}

int Volume_of_Revolution::Ext_Node_ID(int A, int Z)
{
    // Returns the node ID based on the azimuthal and axial position
    if (A>=NA)      return Z*NA + (A-NA);
    else            return Z*NA + A;
}

//----------------------
//------Cylinder--------
//----------------------

void Half_Cylinder::Generate_Nodes()
{
    // We generate the perimeter, and then call the node generation routine from Volume_of_Revolution

    if (Cosine){
        for (int i=0; i<NR; i++)        {Real f = PiCosFac(i,NR+1); Perimeter.push_back(Vector3(R*f,0,Z));        }
        for (int i=0; i<NV+1; i++)      {Real f = PiCosFac(i,NV+1); Perimeter.push_back(Vector3(R,0,Z*(1.0-f)));}
    }
    else{
        for (int i=0; i<NR; i++)    Perimeter.push_back(Vector3(R*i/NR,0,Z));
        for (int i=0; i<NV+1; i++)  Perimeter.push_back(Vector3(R,0,Z*(1.0-1.0*i/NV)));
    }

    Volume_of_Revolution::Generate_Nodes();
}

//-------------------------------
//------Tapered Spar Buoy--------
//-------------------------------

void Tapered_SparBuoy::Generate_Nodes()
{
    // We generate the perimeter, and then call the node generation routine from Volume_of_Revolution
    Real H3 = Z+H1+H2;

    if (Cosine){
        for (int i=0; i<NR-1; i++)      Perimeter.push_back(Vector3(RB*PiCosFac(i,NR),0,Z));
        for (int i=0; i<NV1; i++)       {Real f = PiCosFac(i,NV1+1); Perimeter.push_back(Vector3(RB,0,Z+H1*f));}
        for (int i=0; i<NV2; i++)       {Real f = PiCosFac(i,NV2+1); Perimeter.push_back(Vector3(RB-(RB-RT)*f,0,Z+H1+H2*f));}
        for (int i=0; i<NV3+1; i++)     {Real f = PiCosFac(i,NV3+1); Perimeter.push_back(Vector3(RT,0,H3*(1.0-f)));}
    }
    else{
        for (int i=0; i<NR; i++)        Perimeter.push_back(Vector3(RB*i/NR,0,Z));
        for (int i=0; i<NV1; i++)       Perimeter.push_back(Vector3(RB,0,Z+H1*i/NV1));
        for (int i=0; i<NV2; i++)       Perimeter.push_back(Vector3(RB-(RB-RT)*1.0*i/NV2,0,Z+H1+H2*1.0*i/NV2));
        for (int i=0; i<NV3+1; i++)     Perimeter.push_back(Vector3(RT,0,H3*(1.0-1.0*i/NV3)));
    }

//    for( int i=0; i<Perimeter.size(); i++) std::cout << Perimeter[i](0) <<" "<< Perimeter[i](1) <<" "<< Perimeter[i](2) <<std::endl;

    Volume_of_Revolution::Generate_Nodes();
}

//-------------------------------
//------Sub Leg------------------
//-------------------------------

void Spar_Leg::Generate_Nodes()
{
    // We generate the perimeter, and then call the node generation routine from Volume_of_Revolution
    Real H2 = Z+H1;

    if (Cosine){
        for (int i=0; i<N1; i++)    {Real f = PiCosFac(i,N1+1); Perimeter.push_back(Vector3(RB*f,0,Z));}
        for (int i=0; i<N2; i++)    {Real f = PiCosFac(i,N2+1); Perimeter.push_back(Vector3(RB,0,Z+H1*f));}
        for (int i=0; i<N3; i++)    {Real f = PiCosFac(i,N3+1); Perimeter.push_back(Vector3(RB-(RB-RT)*f,0,Z+H1));}
        for (int i=0; i<N4+1; i++)  {Real f = PiCosFac(i,N4+1); Perimeter.push_back(Vector3(RT,0,H2*(1.0-f)));}
    }
    else{
        for (int i=0; i<N1; i++)    Perimeter.push_back(Vector3(RB*i/N1,0,Z));
        for (int i=0; i<N2; i++)    Perimeter.push_back(Vector3(RB,0,Z+H1*i/N2));
        for (int i=0; i<N3; i++)    Perimeter.push_back(Vector3(RB-(RB-RT)*1.0*i/N3,0,Z+H1));
        for (int i=0; i<N4+1; i++)  Perimeter.push_back(Vector3(RT,0,H2*(1.0-1.0*i/N4)));
    }

    Volume_of_Revolution::Generate_Nodes();
}

}


