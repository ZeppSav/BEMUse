//---------------------------------------------------------------
//----------- Numerical_Wave_Tank Functions----------------------
//---------------------------------------------------------------

#include "Numerical_Wave_Tank.h"

namespace BEMUse
{

//------- Rectangular

void Numerical_Wave_Tank::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Tank_Length"))          L = P.Get_Param<Real>();
        if (P.myNameis("Tank_Width"))           W = P.Get_Param<Real>();
        if (P.myNameis("Tank_Depth"))           D = P.Get_Param<Real>();
        if (P.myNameis("NPanels_Tank_Length"))  NX = P.Get_Param<int>();
        if (P.myNameis("NPanels_Tank_Width"))   NY = P.Get_Param<int>();
        if (P.myNameis("NPanels_Tank_Depth"))   NZ = P.Get_Param<int>();
        if (P.myNameis("Cosine_Disc"))          Cosine = P.Get_Param<bool>();
        if (P.myNameis("Triangular_Panels"))    TriPanels = P.Get_Param<bool>();
    }
}

void Numerical_Wave_Tank::Generate_FreeSurface_Nodes()
{
    // Z Face
    Real x, y;
    for (int i=0; i<NX+1; i++){
        for (int j=0; j<NY+1; j++){
            if (Cosine){
                x = -0.5*L + L*PiCosFac(i,NX+1);
                y = -0.5*W + W*PiCosFac(j,NY+1);
            }
            else{
                x = -0.5*L + L*i/NX;
                y = -0.5*W + W*j/NY;
            }
            Vector3 P(x,y,0);
            Quat tO = Quat::FromTwoVectors(UnitZ,Vector3(0.,0.,-1.));
            FreeSurface_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,P));
        }
    }
}

void Numerical_Wave_Tank::Generate_FreeSurface_Elements()
{
    // Generates new panels
    std::vector<SP_Geo> FSElementtemp;

    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            FSElementtemp.push_back(std::make_shared<Quad_Element>(FreeSurface_Nodes[i*(NY+1)+j],
                                                                   FreeSurface_Nodes[(i+1)*(NY+1)+j],
                                                                   FreeSurface_Nodes[(i+1)*(NY+1)+j+1],
                                                                   FreeSurface_Nodes[i*(NY+1)+j+1]));
        }
    }

    // Convert to triangular panels
    if (TriPanels){
        for (SP_Geo G : FSElementtemp) G->Set_Centroid();    // Ensure Centroids have been created.
        std::vector<SP_Geo> TriPans;
        std::vector<SP_Node> TriNds;
        for (SP_Geo G : FSElementtemp){
            SP_Node C = G->Centroid;
            TriNds.push_back(C);            // Add centroid to array of nodes
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(0),G->Get_Node(1),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(1),G->Get_Node(2),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(2),G->Get_Node(3),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(3),G->Get_Node(0),C));
        }
        StdAppend(FreeSurface_Elements,TriPans);
        StdAppend(FreeSurface_Nodes,TriNds);
    }
    else{
        StdAppend(FreeSurface_Elements,FSElementtemp);
    }

    // Ensure ordering of nodes
    for (int i=0; i<FreeSurface_Nodes.size(); i++) FreeSurface_Nodes[i]->ID = i;
}

void Numerical_Wave_Tank::Generate_Wall_Nodes()
{

    Real x,y,z;

    // X Faces
    for (int j=0; j<NY+1; j++){
        for (int k=0; k<NZ+1; k++){

            if (Cosine){
                y = -0.5*W + W*PiCosFac(j,NY+1);
                z = -D + D*PiCosFac(k,NZ+1);
            }
            else{
                y = -0.5*W + W*j/NY;
                z = -D + D*k/NZ;
            }
            Vector3 P1(-0.5*L,y,z), P2(0.5*L,y,z);
            Quat tb = Quat::FromTwoVectors(UnitZ,Vector3(1.,0.,0.));
            Quat tf = Quat::FromTwoVectors(UnitZ,Vector3(-1.,0.,0.));
            NXFace1.push_back(std::make_shared<Node>(Inertial_CS,tb,P1));
            NXFace2.push_back(std::make_shared<Node>(Inertial_CS,tf,P2));
            // Matrix3 O1 = Matrix3::Zero(); O1(2,0) = -1.; O1(1,1) = 1.; O1(0,2) = 1.;
            // Matrix3 O2 = Matrix3::Zero(); O2(1,0) = -1.; O2(2,1) = 1.; O2(0,2) = -1.;
            // NXFace1.push_back(std::make_shared<Node>(Inertial_CS,Quat(O1),P1));
            // NXFace2.push_back(std::make_shared<Node>(Inertial_CS,Quat(O2),P2));
        }
    }

    // Y Faces
    for (int i=0; i<NX+1; i++){
        for (int k=0; k<NZ+1; k++){

            if (Cosine){
                x = -0.5*L + L*PiCosFac(i,NX+1);
                z = -D + D*PiCosFac(k,NZ+1);
            }
            else{
                x = -0.5*L + L*i/NX;
                z = -D + D*k/NZ;
            }
            Vector3 P1(x,-0.5*W,z), P2(x,0.5*W,z);
            Quat tl = Quat::FromTwoVectors(UnitZ,Vector3(0.,1.,0.));
            Quat tr = Quat::FromTwoVectors(UnitZ,Vector3(0.,-1.,0.));
            NYFace1.push_back(std::make_shared<Node>(Inertial_CS,tl,P1));
            NYFace2.push_back(std::make_shared<Node>(Inertial_CS,tr,P2));
            // Matrix3 O1 = Matrix3::Zero(); O1(2,0) = -1.; O1(0,1) =-1.; O1(1,2) = 1.;
            // Matrix3 O2 = Matrix3::Zero(); O2(2,0) =  1.; O2(0,1) =-1.; O2(1,2) = -1.;
            // NYFace1.push_back(std::make_shared<Node>(Inertial_CS,Quat(O1),P1));
            // NYFace2.push_back(std::make_shared<Node>(Inertial_CS,Quat(O2),P2));
        }
    }

    StdAppend(Wall_Nodes,NXFace1);
    StdAppend(Wall_Nodes,NXFace2);
    StdAppend(Wall_Nodes,NYFace1);
    StdAppend(Wall_Nodes,NYFace2);
}

void Numerical_Wave_Tank::Generate_Wall_Elements()
{
    // Generate wall elements
    std::vector<SP_Geo> WallPans_temp;

    // X Faces
    for (int j=0; j<NY; j++){
        for (int k=0; k<NZ; k++){
            WallPans_temp.push_back(std::make_shared<Quad_Element>(NXFace1[j*(NZ+1)+k],
                                                                   NXFace1[j*(NZ+1)+k+1],
                                                                   NXFace1[(j+1)*(NZ+1)+k+1],
                                                                   NXFace1[(j+1)*(NZ+1)+k]));
            WallPans_temp.push_back(std::make_shared<Quad_Element>(NXFace2[j*(NZ+1)+k],
                                                                   NXFace2[(j+1)*(NZ+1)+k],
                                                                   NXFace2[(j+1)*(NZ+1)+k+1],
                                                                   NXFace2[j*(NZ+1)+k+1]));
        }
    }

    // Y Faces
    for (int i=0; i<NX; i++){
        for (int k=0; k<NZ; k++){
            WallPans_temp.push_back(std::make_shared<Quad_Element>( NYFace1[i*(NZ+1)+k],
                                                                   NYFace1[(i+1)*(NZ+1)+k],
                                                                   NYFace1[(i+1)*(NZ+1)+k+1],
                                                                   NYFace1[i*(NZ+1)+k+1]));
            WallPans_temp.push_back(std::make_shared<Quad_Element>( NYFace2[i*(NZ+1)+k],
                                                                   NYFace2[i*(NZ+1)+k+1],
                                                                   NYFace2[(i+1)*(NZ+1)+k+1],
                                                                   NYFace2[(i+1)*(NZ+1)+k]));
        }
    }

    // Clear node lists (to avoid memory leak)
    NXFace1.clear();
    NXFace2.clear();
    NYFace1.clear();
    NYFace2.clear();

    // Transform panels into triangular panels for RS solver
    if (TriPanels){
        for (SP_Geo G : WallPans_temp) G->Set_Centroid();        // Ensure centroids have been created
        std::vector<SP_Geo> TriPans;
        std::vector<SP_Node> TriNds;
        for (SP_Geo G : WallPans_temp){
            SP_Node C = G->Centroid;
            TriNds.push_back(C);            // Add centroid to array of nodes
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(0),G->Get_Node(1),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(1),G->Get_Node(2),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(2),G->Get_Node(3),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(3),G->Get_Node(0),C));
        }
        StdAppend(Wall_Elements,TriPans);
        StdAppend(Wall_Nodes,TriNds);
    }
    else{
        StdAppend(Wall_Elements,WallPans_temp);
    }

    // Ensure ordering of nodes
    for (int i=0; i<Wall_Nodes.size(); i++) Wall_Nodes[i]->ID = i;
}

void Numerical_Wave_Tank::Generate_Floor_Nodes()
{
    // Z Face
    Real x,y,z;
    for (int i=0; i<NX+1; i++){
        for (int j=0; j<NY+1; j++){
            if (Cosine){
                x = -0.5*L + L*PiCosFac(i,NX+1);     // Cosine
                y = -0.5*W + W*PiCosFac(j,NY+1);     // Cosine
            }
            else{
                x = -0.5*L + L*i/NX;
                y = -0.5*W + W*j/NY;
            }
            Vector3 P(x,y,-D);
            // Floor_Nodes.push_back(std::make_shared<Node>(Inertial_CS,Quat::Identity(),P));
            Quat tO = Quat::FromTwoVectors(UnitZ,Vector3(0.,0.,1.));
            Floor_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,P));
        }
    }
}

void Numerical_Wave_Tank::Generate_Floor_Elements()
{
    // Generate elements of sea bed
    std::vector<SP_Geo> Floor_Elements_temp;

    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            Floor_Elements_temp.push_back(std::make_shared<Quad_Element>(Floor_Nodes[i*(NY+1)+j],
                                                                         Floor_Nodes[i*(NY+1)+j+1],
                                                                         Floor_Nodes[(i+1)*(NY+1)+j+1],
                                                                         Floor_Nodes[(i+1)*(NY+1)+j]));
        }
    }

    if (TriPanels){
        for (SP_Geo G : Floor_Elements_temp) G->Set_Centroid();     // Ensure Centroids have been created.
        std::vector<SP_Geo> TriPans;
        std::vector<SP_Node> TriNds;
        for (SP_Geo G : Floor_Elements_temp){
            SP_Node C = G->Centroid;
            TriNds.push_back(C);            // Add centroid to array of nodes
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(0),G->Get_Node(1),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(1),G->Get_Node(2),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(2),G->Get_Node(3),C));
            TriPans.push_back(std::make_shared<Tri_Element>(G->Get_Node(3),G->Get_Node(0),C));

            // Hack to reorient coordinate system - otherwise we have issues
            Real th = 0.01*D2R;      // Works most of the time(rand()*1.0)/RAND_MAX*M_PI
            int sz = size(TriPans)-1;
            TriPans[sz  ]->Set_Centroid();     TriPans[sz  ]->Centroid->Rotate_Centroid_about_z(th);
            TriPans[sz-1]->Set_Centroid();     TriPans[sz-1]->Centroid->Rotate_Centroid_about_z(th);
            TriPans[sz-2]->Set_Centroid();     TriPans[sz-2]->Centroid->Rotate_Centroid_about_z(th);
            TriPans[sz-3]->Set_Centroid();     TriPans[sz-3]->Centroid->Rotate_Centroid_about_z(th);
        }
        StdAppend(Floor_Elements,TriPans);
        StdAppend(Floor_Nodes,TriNds);
    }
    else{
        StdAppend(Floor_Elements,Floor_Elements_temp);
    }

    // Ensure ordering of nodes
    for (int i=0; i<Floor_Nodes.size(); i++) Floor_Nodes[i]->ID = i;
}

//------- Cylindrical

void Numerical_Wave_Tank_Cylindrical::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Tank_Radius"))              Radius = P.Get_Param<Real>();
        if (P.myNameis("Tank_Depth"))               Depth = P.Get_Param<Real>();
        if (P.myNameis("NPanels_Tank_Azimuth"))     NA = P.Get_Param<int>();
        if (P.myNameis("NPanels_Tank_Radial"))      NR = P.Get_Param<int>();
        if (P.myNameis("NPanels_Tank_Depth"))       NZ = P.Get_Param<int>();
        if (P.myNameis("Cosine_Disc"))              Cosine = P.Get_Param<bool>();
        if (P.myNameis("Triangular_Panels"))        TriPanels = P.Get_Param<bool>();
    }
}

void Numerical_Wave_Tank_Cylindrical::Generate_Wall_Nodes()
{
    // This function generates the walls for the case of a closed geometry
    // The radii of the outside nodes (aligning with the Free-surface & base panels) is taken

    // Values added for convenience.

    if (NZ%2!=0) NZ++;            // Only an EVEN number of vertical sections is allowed
    for (int i=0; i<NZ+1; i++)
    {
        Real Z = -Depth*i/NZ;
        for (int j = 0; j<NA; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NA;
            if (TriPanels) {if (i%2==0) u += 0.5*TwoPI/NA;}   // Shift for triangular panels

            Vector3 P   ( Radius*cos(u), Radius*sin(u), Z);
            Vector3 dx  ( -cos(u), -sin(u), 0.);
            Vector3 dy  ( -sin(u),-cos(u), 0.);
            Vector3 dz  ( 0., 0., 1.);
            Quat tO = Quat::FromTwoVectors(UnitZ,dx);

            Wall_Nodes.push_back(std::make_shared<Node>(Global_CS,tO,P));
        }
    }

    for (int i=0; i<Wall_Nodes.size(); i++) Wall_Nodes[i]->ID = i;
}

void Numerical_Wave_Tank_Cylindrical::Generate_Wall_Elements()
{
    // The wall elements are generated
    // The wall elements are numbered as: [i*NA + j]

    for (int i=0; i<NZ; i++)
    {
        for (int j = 0; j<NA; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            int a1 = j;
            int a2 = j+1;
            if (a2>=NA) a2 -= NA;

            if (i%2==0){
                Wall_Elements.push_back(std::make_shared<Tri_Element>(Wall_Nodes[i*NA + a1],
                                                                      Wall_Nodes[(i+1)*NA + a2],
                                                                      Wall_Nodes[i*NA + a2]));
                Wall_Elements.push_back(std::make_shared<Tri_Element>(Wall_Nodes[i*NA + a1],
                                                                      Wall_Nodes[(i+1)*NA + a1],
                                                                      Wall_Nodes[(i+1)*NA + a2]));
            }
            else
            {
                Wall_Elements.push_back(std::make_shared<Tri_Element>(Wall_Nodes[i*NA + a1],
                                                                      Wall_Nodes[(i+1)*NA + a1],
                                                                      Wall_Nodes[i*NA + a2]));

                Wall_Elements.push_back(std::make_shared<Tri_Element>(Wall_Nodes[i*NA + a2],
                                                                      Wall_Nodes[(i+1)*NA + a1],
                                                                      Wall_Nodes[(i+1)*NA + a2]));
            }
        }
    }

    for (int i=0; i<Wall_Elements.size(); i++)  Wall_Elements[i]->Set_Centroid();
}

void Numerical_Wave_Tank_Cylindrical::Generate_Floor_Nodes()
{
    // We generate here the node for the Free-surface grid. using triangular panels

    for (int i=0; i<NR+1; i++)
    {
        Real r;
        if (Cosine) r = Radius*(1.0-HalfPiCosFac(NR-i,NR+1));
        else        r = Radius*1.0*i/NZ;

        int NAR = NA;
        if (i==0)      NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;

            //--- Mod for triangular panels
            if (TriPanels) if (i%2==0) u += 0.5*TwoPI/NAR;

            Vector3 P   ( r*cos(u), r*sin(u), -Depth);
            Vector3 dx  ( sin(u), cos(u), 0.);
            Vector3 dy  ( -cos(u),sin(u), 0.);
            Vector3 dz  ( 0., 0., 1.);
            Quat tO = Quat::FromTwoVectors(UnitX,dx);
            // Quat tO = Quat::FromTwoVectors(UnitZ,-dz);

            Floor_Nodes.push_back(std::make_shared<Node>(Global_CS,tO,P));
        }
    }

    for (int i=0; i<Floor_Nodes.size(); i++) Floor_Nodes[i]->ID = i;
}

void Numerical_Wave_Tank_Cylindrical::Generate_Floor_Elements()
{
    // This function generates external (free-surface) panels

    // std::cout << "Generating Floor Elements \n";

    // Base (triangular elements)
    for (int i=0; i<NA; i++)    Floor_Elements.push_back(std::make_shared<Tri_Element>( Floor_Nodes[Ext_Node_ID(0,0)],
                                                                                        Floor_Nodes[Ext_Node_ID(i+1,1)],
                                                                                        Floor_Nodes[Ext_Node_ID(i,1)]));

    // Generate central elements with triangular panels
    for (int z=1; z<NR; z++){
        for (int i=0; i<NA; i++){
            if (z%2==0){
                Floor_Elements.push_back(std::make_shared<Tri_Element>( Floor_Nodes[Ext_Node_ID(i,z)],
                                                                       Floor_Nodes[Ext_Node_ID(i+1,z)],
                                                                       Floor_Nodes[Ext_Node_ID(i+1,z+1)]));
                Floor_Elements.push_back(std::make_shared<Tri_Element>( Floor_Nodes[Ext_Node_ID(i+1,z)],
                                                                       Floor_Nodes[Ext_Node_ID(i+2,z+1)],
                                                                       Floor_Nodes[Ext_Node_ID(i+1,z+1)]));
            }
            else
            {
                Floor_Elements.push_back(std::make_shared<Tri_Element>( Floor_Nodes[Ext_Node_ID(i,z)],
                                                                       Floor_Nodes[Ext_Node_ID(i+1,z)],
                                                                       Floor_Nodes[Ext_Node_ID(i,z+1)]));
                Floor_Elements.push_back(std::make_shared<Tri_Element>( Floor_Nodes[Ext_Node_ID(i+1,z)],
                                                                       Floor_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                       Floor_Nodes[Ext_Node_ID(i,z+1)]));
            }
        }
    }

    for (int i=0; i<Floor_Elements.size(); i++)  Floor_Elements[i]->Set_Centroid();
}

int Numerical_Wave_Tank_Cylindrical::Ext_Node_ID(int A, int Z)
{
    // Returns the id of external (free surface) grid nodes
    if (Z==0)      return 0;
    else {
        if (A>=NA)    return 1+(Z-1)*NA+(A-NA);
        else            return 1+(Z-1)*NA+A;
    }
}

void Numerical_Wave_Tank_Cylindrical::Generate_FreeSurface_Nodes()
{
    // We generate here the node for the Free-surface grid. using triangular panels

    for (int i=0; i<NR+1; i++)
    {
        Real r;
        if (Cosine) r = Radius*(1.0-HalfPiCosFac(NR-i,NR+1));
        else        r = Radius*1.0*i/NZ;

        int NAR = NA;
        if (i==0)      NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;

            //--- Mod for triangular panels
            if (TriPanels) {if (i%2==0) u += 0.5*TwoPI/NAR;}

            Vector3 P   ( r*cos(u), r*sin(u), 0.);
            Vector3 dx  ( sin(u), cos(u), 0.);
            Vector3 dy  ( -cos(u),sin(u), 0.);
            Vector3 dz  ( 0., 0., 1.);
            // Quat tO = Quat::FromTwoVectors(UnitX,dx);
            Quat tO = Quat::FromTwoVectors(UnitZ,-dz);

            FreeSurface_Nodes.push_back(std::make_shared<Node>(Global_CS,tO,P));
        }
    }

    for (int i=0; i<FreeSurface_Nodes.size(); i++) FreeSurface_Nodes[i]->ID = i;
}

void Numerical_Wave_Tank_Cylindrical::Generate_FreeSurface_Elements()
{
    // This function generates external (free-surface) panels

    // std::cout << "Generating Free-Surface Elements \n";

    // Base (triangular elements)
    for (int i=0; i<NA; i++)    FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(   FreeSurface_Nodes[Ext_Node_ID(0,0)],
                                                                                                FreeSurface_Nodes[Ext_Node_ID(i,1)],
                                                                                                FreeSurface_Nodes[Ext_Node_ID(i+1,1)]));

    // Generate central elements with triangular panels
    for (int z=1; z<NR; z++){
        for (int i=0; i<NA; i++){
            if (z%2==0){
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i,z)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i+1,z)]));
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i+1,z)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i+2,z+1)]));
            }
            else
            {
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i,z)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i,z+1)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i+1,z)]));
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i+1,z)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i,z+1)],
                                                                             FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)]));
            }
        }
    }

    for (int i=0; i<FreeSurface_Elements.size(); i++)  FreeSurface_Elements[i]->Set_Centroid();
}


}
