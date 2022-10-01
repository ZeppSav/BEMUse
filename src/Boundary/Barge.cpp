//---------------------------------------------------------------
//------------------ Barge Functions-----------------------------
//---------------------------------------------------------------

#include "Barge.h"

namespace BEMUse
{

//--- Geometry functions

void Barge::Generate_Nodes()
{
    // For simplicity I do this over 5 faces. These are then discarded in the element creation module to avoid issues upon deletion.

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

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
            NXFace1.push_back(std::make_shared<Node>(Inertial_CS,P1));
            NXFace2.push_back(std::make_shared<Node>(Inertial_CS,P2));
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
            NYFace1.push_back(std::make_shared<Node>(Inertial_CS,P1));
            NYFace2.push_back(std::make_shared<Node>(Inertial_CS,P2));
        }
    }

    // Z Face
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
            NZFace.push_back(std::make_shared<Node>(Inertial_CS,P));
        }
    }

    StdAppend(Nodes,NXFace1);
    StdAppend(Nodes,NXFace2);
    StdAppend(Nodes,NYFace1);
    StdAppend(Nodes,NYFace2);
    StdAppend(Nodes,NZFace);

    for (int i=0; i<Nodes.size(); i++) Nodes[i]->ID = i;
}

void Barge::Generate_Aux_Nodes()
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
                y = -0.5*W + W*i/NY;
            }
            Vector3 P(x,y,0);
            Aux_Nodes.push_back(std::make_shared<Node>(Inertial_CS,P));
        }
    }

    for (int i=0; i<Aux_Nodes.size(); i++) Aux_Nodes[i]->ID = i;
}

void Barge::Generate_Ext_Nodes()
{
    Real x, y;
    // For simplicity I will simply generate a single square grid, and later I will remove the unnecessary elements
    for (int i=0; i<NXFS+1; i++){
        for (int j=0; j<NYFS+1; j++){
            x = -0.5*LFS + LFS*i/NXFS;
            y = -0.5*WFS + WFS*j/NYFS;      // Only linear here!
            // As a trick to extend the nodes SLIGHTLY away from the boundary I will dilate them here
            Vector3 P(1.001*x,1.001*y,0);
            Ext_Nodes.push_back(std::make_shared<Node>(Inertial_CS,P));
        }
    }

    for (int i=0; i<Ext_Nodes.size(); i++) Ext_Nodes[i]->ID = i;
}

void Barge::Generate_Elements()
{
    // All geo elements will be define here using the Previously specified nodes.

    // X Faces
    for (int j=0; j<NY; j++){
        for (int k=0; k<NZ; k++){
            Elements.push_back(std::make_shared<Quad_Element>(  NXFace1[j*(NZ+1)+k],
                                                                NXFace1[j*(NZ+1)+k+1],
                                                                NXFace1[(j+1)*(NZ+1)+k+1],
                                                                NXFace1[(j+1)*(NZ+1)+k]));
            Elements.push_back(std::make_shared<Quad_Element>(  NXFace2[j*(NZ+1)+k],
                                                                NXFace2[(j+1)*(NZ+1)+k],
                                                                NXFace2[(j+1)*(NZ+1)+k+1],
                                                                NXFace2[j*(NZ+1)+k+1]));
        }
    }

    // Y Faces
    for (int i=0; i<NX; i++){
        for (int k=0; k<NZ; k++){
            Elements.push_back(std::make_shared<Quad_Element>(  NYFace1[i*(NZ+1)+k],
                                                                NYFace1[(i+1)*(NZ+1)+k],
                                                                NYFace1[(i+1)*(NZ+1)+k+1],
                                                                NYFace1[i*(NZ+1)+k+1]));
            Elements.push_back(std::make_shared<Quad_Element>(  NYFace2[i*(NZ+1)+k],
                                                                NYFace2[i*(NZ+1)+k+1],
                                                                NYFace2[(i+1)*(NZ+1)+k+1],
                                                                NYFace2[(i+1)*(NZ+1)+k]));
        }
    }

    // Z Face
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            Elements.push_back(std::make_shared<Quad_Element>(  NZFace[i*(NY+1)+j],
                                                                NZFace[i*(NY+1)+j+1],
                                                                NZFace[(i+1)*(NY+1)+j+1],
                                                                NZFace[(i+1)*(NY+1)+j]));
        }
    }

    // Clear node lists (to avoid memory leak)
    NXFace1.clear();
    NXFace2.clear();
    NYFace1.clear();
    NYFace2.clear();
    NZFace.clear();

}

void Barge::Generate_Aux_Elements()
{
    // Generates free surface panels (Irregular frequency removal)
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            Aux_Elements.push_back(std::make_shared<Quad_Element>(  Aux_Nodes[i*(NY+1)+j],
                                                                    Aux_Nodes[(i+1)*(NY+1)+j],
                                                                    Aux_Nodes[(i+1)*(NY+1)+j+1],
                                                                    Aux_Nodes[i*(NY+1)+j+1]));
        }
    }
}

void Barge::Generate_Ext_Elements()
{
    // Create elements, everywhere except WITHIN the interior free surface

    std::vector<bool> Applied;
    Applied.assign(Ext_Nodes.size(),false);

    Real XL = -0.5*L, XR = 0.5*L;
    Real YL = -0.5*W, YR = 0.5*W;

    for (int i=0; i<NXFS; i++){
        for (int j=0; j<NYFS; j++){
            SP_Node N1 = Ext_Nodes[i*(NYFS+1)+j];           Vector3 P1 = N1->Position_Global();
            SP_Node N2 = Ext_Nodes[(i+1)*(NYFS+1)+j];       Vector3 P2 = N2->Position_Global();
            SP_Node N3 = Ext_Nodes[(i+1)*(NYFS+1)+j+1];     Vector3 P3 = N3->Position_Global();
            SP_Node N4 = Ext_Nodes[i*(NYFS+1)+j+1];         Vector3 P4 = N4->Position_Global();

            Vector3 PMid = 0.25*(P1+P2+P3+P4);
            bool InX = (fabs(PMid(0))<0.5*L);
            bool InY = (fabs(PMid(1))<0.5*W);
            if (InX&&InY) continue;             // Do not include panels over the interior free surface

            Ext_Elements.push_back(std::make_shared<Quad_Element>(N1,N2,N3,N4));

            // Set these nodes to be active
            Applied[N1->ID] = true;
            Applied[N2->ID] = true;
            Applied[N3->ID] = true;
            Applied[N4->ID] = true;

        }
    }

    // Create surrogate list to store additional values
    std::vector<SP_Node> Surr_Nodes;
    for (int i=0; i<Ext_Nodes.size(); i++) if (Applied[i]) Surr_Nodes.push_back(Ext_Nodes[i]);

    // Now replace array (& reset node id)
    Ext_Nodes.clear();
    StdAppend(Ext_Nodes,Surr_Nodes);
    for (int i=0; i<Ext_Nodes.size(); i++) Ext_Nodes[i]->ID = i;
}

}
