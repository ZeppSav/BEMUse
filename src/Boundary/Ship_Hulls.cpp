//---------------------------------------------------------------
//------------------ Ship Hull Functions-------------------------
//---------------------------------------------------------------

#include "Ship_Hulls.h"

namespace BEMUse
{

//--- Generate elements
void Ship_Hull::Generate_Elements()
{
    // Creates panels for the wigley hull

    std::vector<SP_Geo> LeftPans, RightPans;

    for (int i=0; i<NX; i++)
    {
        for (int j=0; j<NZ; j++)
        {
            LeftPans.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,j,0)],
                                                                Nodes[Node_ID(i,j+1,0)],
                                                                Nodes[Node_ID(i+1,j+1,0)],
                                                                Nodes[Node_ID(i+1,j,0)]));
            RightPans.push_back(std::make_shared<Quad_Element>( Nodes[Node_ID(i+1,j,1)],
                                                                Nodes[Node_ID(i+1,j+1,1)],
                                                                Nodes[Node_ID(i,j+1,1)],
                                                                Nodes[Node_ID(i,j,1)]));
        }
    }

    StdAppend(Elements,LeftPans);
    StdAppend(Elements,RightPans);

    LeftPans.clear();
    RightPans.clear();
}

//--- ID retrieval
int Ship_Hull::Node_ID(int X, int Z, int L)
{
    // Elements of an area
    return L*(NX+1)*(NZ+1) + X*(NZ+1)+Z;
}

//----------------------
//------Wigley Hull-----
//----------------------

void Wigley_Hull::Generate_Nodes()
{
    // Creates Nodes for the wigley hull

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    std::vector<SP_Node> LeftNodes, RightNodes;
    Real x,z;

    for (int i=0; i<NX+1; i++)
    {
        for (int j=0; j<NZ+1; j++)
        {
            if (Cosine){
                x = L*PiCosFac(i,NX+1);
                z = T*(PiCosFac(j,NZ+1)-1.0);
            }
            else{
                x = L*i/NX;
                z = -T + T*j/NZ;
            }

            Real xf = x/L, zf = z/T;
            Real y = 2.0*B*(xf*(1.0-xf))*(1.0-zf*zf);

            Vector3 BPL(x-0.5*L,y,z), BPR(x-0.5*L,-y,z);
            LeftNodes.push_back(std::make_shared<Node>(Inertial_CS,BPL));
            RightNodes.push_back(std::make_shared<Node>(Inertial_CS,BPR));
        }
    }

    StdAppend(Nodes,LeftNodes);
    StdAppend(Nodes,RightNodes);

    LeftNodes.clear();
    RightNodes.clear();

    // Set node ID
    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;
}

}
