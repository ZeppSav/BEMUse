//---------------------------------------------------------------
//------------------ Thin Disc Functions-------------------------
//---------------------------------------------------------------

#include "Thin_Disc.h"

namespace BEMUse
{

void Thin_Disc::Generate_Nodes()
{
    // Simplest possible geo. A thin disc

    if (!Global_CS)    Global_CS = new CoordSys();
//    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    Vector3 P(0,0,Depth);
    Quat O(Quat::Identity());
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS,O,P); // Rotated

    // Create nodes
    Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(0.0,0.0,0.0)));      // Centre node

    for (int i=1; i<NR+1; i++){
        Real r;
        if (Cosine) r = R*PiCosFac(i,NR+1);
        else        r = R*i/(NR-1);
        for (int j=0; j<NA; j++){
            Real th = j*TwoPI/NA;
            Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(r*sin(th),r*cos(th),0.0)));
        }
    }
}

void Thin_Disc::Generate_Elements()
{
    // Centre panels
    for (int j=0; j<NA; j++)    Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                                                    Nodes[Node_ID(j,1)],
                                                                                    Nodes[Node_ID(j+1,1)]));

    // Outer panels
    for (int i=1; i<NR; i++){
         for (int j=0; j<NA; j++) Elements.push_back(std::make_shared<Quad_Element>(Nodes[Node_ID(j,i)],
                                                                                    Nodes[Node_ID(j,i+1)],
                                                                                    Nodes[Node_ID(j+1,i+1)],
                                                                                    Nodes[Node_ID(j+1,i)] ));
    }

//    for (SP_Geo G : Elements)    G->isDipole = true;
}

int Thin_Disc::Node_ID(int A, int Z)
{
    if (Z==0) return 0;
    else {
        if (A>=NA)      return 1+(Z-1)*NA+(A-NA);
        else if (A<0)   return 1+(Z-1)*NA+(A+NA);
        else            return 1+(Z-1)*NA+A;
    }
}

}
