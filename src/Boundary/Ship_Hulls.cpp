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

    if (TriPanels)  Generate_Triangular_Elements();
    else            Generate_Quadratic_Elements();

    // Ensure ordering of nodes
    // for (size_t i=0; i<size(Nodes); i++) Nodes[i]->ID = i;
}

void Ship_Hull::Generate_Quadratic_Elements()
{
    // Creates panels for the wigley hull

    std::vector<SP_Geo> LeftPans, RightPans;

    for (int i=0; i<NX; i++)
    {
        for (int j=0; j<NZ; j++)
        {
            LeftPans.push_back(std::make_shared<Quad_Element>(  LeftNodes[ i*(NZ+1)     +   j],
                                                              LeftNodes[ i*(NZ+1)     +   j+1],
                                                              LeftNodes[ (i+1)*(NZ+1) +   j+1],
                                                              LeftNodes[ (i+1)*(NZ+1) +   j]));

            RightPans.push_back(std::make_shared<Quad_Element>( RightNodes[(i+1)*(NZ+1) +   j],
                                                               RightNodes[(i+1)*(NZ+1) +   j+1],
                                                               RightNodes[i*(NZ+1)     +   j+1],
                                                               RightNodes[i*(NZ+1)     +   j]));
        }
    }

    StdAppend(Elements,LeftPans);
    StdAppend(Elements,RightPans);

    // Reverse ordering if we wish for outward pointing normals
    if (Panels_Outward){
        for (SP_Geo G : Elements) G->Reorder_Nodes();
    }
}

void Ship_Hull::Generate_Triangular_Elements()
{
    // Creates panels for the wigley hull

    std::vector<SP_Geo> LeftPans, RightPans;

    for (int i=0; i<NX; i++)
    {
        for (int j=0; j<NZ; j++)
        {
            // Left panel
            SP_Node LN1 = LeftNodes[ i*(NZ+1)     +   j     ];
            SP_Node LN2 = LeftNodes[ i*(NZ+1)     +   j+1   ];
            SP_Node LN3 = LeftNodes[ (i+1)*(NZ+1) +   j+1   ];
            SP_Node LN4 = LeftNodes[ (i+1)*(NZ+1) +   j     ];
            SP_Node LC = TriLeftNodes[i*NZ + j];

            LeftPans.push_back(std::make_shared<Tri_Element>(LN1, LN2, LC));
            LeftPans.push_back(std::make_shared<Tri_Element>(LN2, LN3, LC));
            LeftPans.push_back(std::make_shared<Tri_Element>(LN3, LN4, LC));
            LeftPans.push_back(std::make_shared<Tri_Element>(LN4, LN1, LC));

            SP_Node RN1 = RightNodes[(i+1)*(NZ+1) +   j];
            SP_Node RN2 = RightNodes[(i+1)*(NZ+1) +   j+1];
            SP_Node RN3 = RightNodes[i*(NZ+1)     +   j+1];
            SP_Node RN4 = RightNodes[i*(NZ+1)     +   j];
            SP_Node RC = TriRightNodes[i*NZ + j];

            RightPans.push_back(std::make_shared<Tri_Element>(RN1, RN2, RC));
            RightPans.push_back(std::make_shared<Tri_Element>(RN2, RN3, RC));
            RightPans.push_back(std::make_shared<Tri_Element>(RN3, RN4, RC));
            RightPans.push_back(std::make_shared<Tri_Element>(RN4, RN1, RC));
        }
    }

    StdAppend(Elements,LeftPans);
    StdAppend(Elements,RightPans);

    // Reverse ordering if we wish for outward pointing normals
    if (Panels_Outward){
        for (SP_Geo G : Elements) G->Reorder_Nodes();
    }
}

//----------------------
//------Wigley Hull-----
//----------------------

void Wigley_Hull::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Length"))               L = P.Get_Param<Real>();
        if (P.myNameis("Beam"))                 B = P.Get_Param<Real>();
        if (P.myNameis("Draft"))                T = P.Get_Param<Real>();
        if (P.myNameis("NPanels_Length"))       NX = P.Get_Param<int>();
        if (P.myNameis("NPanels_Depth"))        NZ = P.Get_Param<int>();
        if (P.myNameis("NPanels_FS_Length_Bow"))    NXS = P.Get_Param<int>();
        if (P.myNameis("NPanels_FS_Length_Stern"))  NXK = P.Get_Param<int>();
        if (P.myNameis("NPanels_FS_Width"))     NYFS = P.Get_Param<int>();
        if (P.myNameis("Cosine_Disc"))          Cosine = P.Get_Param<bool>();
        if (P.myNameis("Triangular_Panels"))    TriPanels = P.Get_Param<bool>();
        if (P.myNameis("Panels_Outward"))       Panels_Outward = P.Get_Param<bool>();

        if (P.myNameis("BowGridFactor"))        XFSMax = P.Get_Param<Real>()*L;
        if (P.myNameis("SternGridFactor"))      XFSMin = -P.Get_Param<Real>()*L;
        if (P.myNameis("PortGridFactor"))       YFS = P.Get_Param<Real>()*L;
    }
    // TriPanels = true;       // HACK!

    // Set default parameters in case not already set
    XHK = -0.5*L;
    XHS = 0.5*L;
    if (NXS == 0)       NXS = NX;
    if (NXK == 0)       NXK = 2*NX;
    if (NYFS == 0)      NYFS = 1*NZ;
    if (XFSMax == 0)    XFSMax = 1.0*L;
    if (XFSMin == 0)    XFSMin = -2.0*L;
    if (YFS == 0)       YFS = L;
}

void Wigley_Hull::Hull_Node(Real xn, Real zn, Vector3 &p, Vector3 &n)
{
    // This function takes in the normlized hull coordinates xn in [0, 1], zn in [-1,0] and returns the node position

    p(0) = L*xn - 0.5*L;                    // Shift to ensure centre beam at origin
    p(1) = 2.0*B*xn*(1.0-xn)*(1.0-zn*zn);
    p(2) = T*zn;

    // Set normal
    // Hull surface given by:
    // f = y - 2B xf (1-xf) (1-zf^2)
    // Gradient given by:
    // df/dx = -2B (1-zf^2)1/L (1 - 2x/L)
    // df/dy = 1
    // df/dz = 4B xf (1-xf) 2 zf/t
    n(0) = -2.0*B*(1.0-zn*zn)*(1.0-2.0*xn)/L;
    n(1) = 1.0;
    n(2) = 4.0*B*xn*(1.0-xn)*zn/T;
    n.normalize();
}

void Wigley_Hull::Generate_Nodes()
{
    // Creates Nodes for the wigley hull

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    // std::vector<SP_Node> LeftNodes, RightNodes;
    Real xf,zf;
    Vector3 pn, nn;

    for (int i=0; i<NX+1; i++)
    {
        if (Cosine) xf = 1.0*PiCosFac(i,NX+1);
        else        xf = 1.0*i/NX;

        Hull_Node(xf,0,pn,nn);
        Hull_Perim.push_back(pn);

        for (int k=0; k<NZ+1; k++)
        {

            if (Cosine) zf = PiCosFac(k,NZ+1)-1.0;
            else        zf = -1. + 1.0*k/NZ;

            // std::cout << "Main " << zf << std::endl;

            Hull_Node(xf,zf,pn,nn);
            Quat tOL = Quat::FromTwoVectors(UnitZ,nn);
            LeftNodes.push_back(std::make_shared<Node>(Inertial_CS, tOL, pn));

            // Switch position/normal for node on opposite side
            pn(1) *= -1;
            nn(1) *= -1;
            Quat tOR = Quat::FromTwoVectors(UnitZ,nn);
            RightNodes.push_back(std::make_shared<Node>(Inertial_CS, tOR, pn));
        }
    }

    StdAppend(Nodes,LeftNodes);
    StdAppend(Nodes,RightNodes);

    // Generate panel centre nodes in the case there are triangular panels
    // This approach is a little cleaner than using the centroid of the generated panels...

    if (TriPanels){

        for (int i=0; i<NX; i++)
        {
            if (Cosine) xf = 1.0*PiCosFacShift(i,NX+1);
            else        xf = (i+0.5)/NX;

            for (int k=0; k<NZ; k++)
            {

                if (Cosine) zf = PiCosFacShift(k,NZ+1)-1.0;
                else        zf = -1. + (k+0.5)/NZ;

                Hull_Node(xf,zf,pn,nn);
                Quat tOL = Quat::FromTwoVectors(UnitZ,nn);
                TriLeftNodes.push_back(std::make_shared<Node>(Inertial_CS, tOL, pn));

                // Switch for node on opposite position
                pn(1) *= -1;
                nn(1) *= -1;
                Quat tOR = Quat::FromTwoVectors(UnitZ,nn);
                TriRightNodes.push_back(std::make_shared<Node>(Inertial_CS, tOR, pn));
            }
        }

        StdAppend(Nodes,TriLeftNodes);
        StdAppend(Nodes,TriRightNodes);

    }

    // Set node ID
    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;
    // std::cout << "Body Nodes" << std::endl;
    // for (int i=0; i<size(Nodes); i++) std::cout <<     Nodes[i]->ID << std::endl;
}

void Wigley_Hull::Generate_Aux_Nodes()
{
    // Creates auxiliary nodes for the interior free surface of the wigley hull
    // I shall assume for now that the same discretisation is being used as the ship body

    for (int i=0; i<Hull_Perim.size(); i++)
    {
        for (int j=0; j<NZ+1; j++)
        {
            Real f;
            if (Cosine) f = HalfPiCosFac(j,NZ+1);
            else        f = 1.0*j/NZ;

            Real x = Hull_Perim[i](0);
            Real y = f*Hull_Perim[i](1);
            Aux_Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(x,y,0.0)));
        }
    }

    // Set node ID
    for (int i=0; i<Aux_Nodes.size(); i++)  Aux_Nodes[i]->ID = i;
}

void Wigley_Hull::FS_Node(Real xn, Real yn, Vector3 &p, int zone)
{
    // This function takes in the normlized FS coordinates xn in [0, 1], zn in [0, 1] and returns the node position
    // Zone represents (with ship translating in positive x-direction) 0: Aft, 1: beam, 2: Fore

    // What are the dimensions of the FS grid?
    NXH = Hull_Perim.size()-1;

    if (zone==0){       // FS Region aft of boat
        p(0) = XFSMin + (XHK-XFSMin)*xn;
        p(1) = YFS*yn;
        p(2) = 0;
    }

    if (zone==1){       // FS Region next to of boat
        Real ymin = 2.0*B*xn*(1.0-xn);      // y-coord at the top of the wigley hull
        p(0) =  -0.5*L + L*xn;
        p(1) = ymin + (YFS-ymin)*yn;
        p(2) = 0;
    }

    if (zone==2){       // FS Region fore of boat
        p(0) = XHS + (XFSMax-XHS)*xn;
        p(1) = YFS*yn;
        p(2) = 0;
    }
}

void Wigley_Hull::Generate_FreeSurface_Nodes()
{
    // Using the Hull_Perim variable we will create the X grid.
    // First create a list of X/ Y positions for the starting points

    // Specify dummy values
    NXH = Hull_Perim.size()-1;
    XHK = -0.5*L;
    XHS = 0.5*L;

    // The non-dimensional (x,y) coordinates are stored first.

    Real fx, fy;
    Vector3 pn;
    std::vector<Real> FSXGrid1, FSXGrid2, FSXGrid3, FSYGrid1, FSYGrid2, FSYGrid3;

    // X-arrays

    // Generate section aft of boat
    for (int i=0; i<NXK; i++){
        if (Cosine) fx = HalfPiCosFac(i,NXK+1);
        else        fx = 1.0*i/NXK;
        FSXGrid1.push_back(fx);
    }
    // Generate section beam of boat
    for (int i=0; i<NX+1; i++){
        if (Cosine) fx = PiCosFac(i,NX+1);
        else        fx = 1.0*i/NX;
        FSXGrid2.push_back(fx);
    }
    // Generate section fore of boat
    for (int i=1; i<NXS+1; i++){
        if (Cosine) fx = HalfPiCosFac2(i,NXS+1);      // HERE
        else        fx = 1.0*i/NXS;
        FSXGrid3.push_back(fx);
    }

    // Y-arrays
    for (int j=0; j<NYFS+1; j++){
        if (Cosine) fy = HalfPiCosFac2(j,NYFS+1);
        else        fy = 1.0*j/NYFS;
        FSYGrid1.push_back(fy);
        // FSYGrid2.push_back(fy);
        // Hack to avoid singularities at hull intersection
        if (j==0)   FSYGrid2.push_back(0.001);
        else        FSYGrid2.push_back(fy);
        FSYGrid3.push_back(fy);
    }

    // Now generate nodes
    Quat tO = Quat::FromTwoVectors(UnitZ,Vector3(0,0,-1.));
    for (size_t i=0; i<size(FSXGrid1); i++){
        for (int j=0; j<NYFS+1; j++){
            FS_Node(FSXGrid1[i],FSYGrid1[j],pn,0);
            FreeSurface_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,pn));
        }
    }
    for (size_t i=0; i<size(FSXGrid2); i++){
        for (int j=0; j<NYFS+1; j++){
            FS_Node(FSXGrid2[i],FSYGrid2[j],pn,1);
            FreeSurface_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,pn));
        }
    }
    for (size_t i=0; i<size(FSXGrid3); i++){
        for (int j=0; j<NYFS+1; j++){
            FS_Node(FSXGrid3[i],FSYGrid3[j],pn,2);
            FreeSurface_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,pn));
        }
    }

    // for (int i=0; i<size(FreeSurface_Nodes); i++){
    //     Vector3 G = FreeSurface_Nodes[i]->Position_Global();
    //     std::cout <<  G(0) csp G(1) csp G(2)  << std::endl;
    // }

    if (TriPanels){

        // Clear arrays from above
        FSXGrid1.clear();
        FSXGrid2.clear();
        FSXGrid3.clear();
        FSYGrid1.clear();
        FSYGrid2.clear();
        FSYGrid3.clear();

        // Generate section aft of boat
        for (int i=0; i<NXK; i++){
            if (Cosine) fx = HalfPiCosFacShift(i,NXK+1);
            else        fx = 1.0*(i+0.5)/NXK;
            FSXGrid1.push_back(fx);
        }
        // Generate section beam of boat
        for (int i=0; i<NX; i++){
            if (Cosine) fx = PiCosFacShift(i,NX+1);
            else        fx = 1.0*(i+0.5)/NX;
            FSXGrid2.push_back(fx);
        }
        // Generate section fore of boat
        for (int i=0; i<NXS; i++){
            if (Cosine) fx = HalfPiCosFac2Shift(i,NXS+1);
            else        fx = 1.0*(i+0.5)/(NXS);
            FSXGrid3.push_back(fx);
        }

        // Y-arrays

        for (int j=0; j<NYFS; j++){
            if (Cosine) fy = HalfPiCosFac2Shift(j,NYFS+1);
            else        fy = 1.0*(j+0.5)/NYFS;
            FSYGrid1.push_back(fy);
            FSYGrid2.push_back(fy);
            FSYGrid3.push_back(fy);
        }

        // Now generate nodes
        Quat tO = Quat::FromTwoVectors(UnitZ,Vector3(0,0,-1.));
        for (size_t i=0; i<size(FSXGrid1); i++){
            for (int j=0; j<NYFS; j++){
                FS_Node(FSXGrid1[i],FSYGrid1[j],pn,0);
                Tri_FS_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,pn));
            }
        }
        for (size_t i=0; i<size(FSXGrid2); i++){
            for (int j=0; j<NYFS; j++){
                FS_Node(FSXGrid2[i],FSYGrid2[j],pn,1);
                Tri_FS_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,pn));
            }
        }
        for (size_t i=0; i<size(FSXGrid3); i++){
            for (int j=0; j<NYFS; j++){
                FS_Node(FSXGrid3[i],FSYGrid3[j],pn,2);
                Tri_FS_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,pn));
            }
        }

        StdAppend(FreeSurface_Nodes,Tri_FS_Nodes);

        // for (int i=0; i<size(Tri_FS_Nodes); i++){
        //     Vector3 G = Tri_FS_Nodes[i]->Position_Global();
        //     std::cout <<  G(0) csp G(1) csp G(2)  << std::endl;
        // }
    }
}

void Wigley_Hull::Generate_Aux_Elements()
{
    // Creates panels for the wigley hull

    for (int i=0; i<Hull_Perim.size()-1; i++){
        for (int j=0; j<NZ; j++)
        {
            Aux_Elements.push_back(std::make_shared<Quad_Element>(  Aux_Nodes[i*(NZ+1) + j],
                                                                  Aux_Nodes[i*(NZ+1) + j+1],
                                                                  Aux_Nodes[(i+1)*(NZ+1) + j+1],
                                                                  Aux_Nodes[(i+1)*(NZ+1) + j]));
        }
    }
}

void Wigley_Hull::Generate_FreeSurface_Elements()
{
    // Creates panels for the wigley hull

    int NXT = NXK + NXH + NXS;
    // int NXT = NXK + NXH + NXS - 2;

    if (TriPanels){         // Triangular panels

        for (int i=0; i<NXT; i++){
            for (int j=0; j<NYFS; j++)
            {
                SP_Node N1 = FreeSurface_Nodes[i*(NYFS+1) + j];
                SP_Node N2 = FreeSurface_Nodes[i*(NYFS+1) + j+1];
                SP_Node N3 = FreeSurface_Nodes[(i+1)*(NYFS+1) + j+1];
                SP_Node N4 = FreeSurface_Nodes[(i+1)*(NYFS+1) + j];
                SP_Node C = Tri_FS_Nodes[i*NYFS + j];
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(N1, N2, C));
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(N2, N3, C));
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(N3, N4, C));
                FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(N4, N1, C));
            }
        }
    }
    else {                  // Quadrilateral panels

        for (int i=0; i<NXT; i++){
            for (int j=0; j<NYFS; j++)
            {
                FreeSurface_Elements.push_back(std::make_shared<Quad_Element>(  FreeSurface_Nodes[i*(NYFS+1) + j],
                                                                              FreeSurface_Nodes[i*(NYFS+1) + j+1],
                                                                              FreeSurface_Nodes[(i+1)*(NYFS+1) + j+1],
                                                                              FreeSurface_Nodes[(i+1)*(NYFS+1) + j]));
            }
        }

    }

    // std::cout << "N Free surface nodes " << size(FreeSurface_Nodes) csp NXT+1 csp NYFS+1 << std::endl;
    // std::cout << "N Free surface panels " << size(FreeSurface_Elements) << std::endl;

    // Ensure ordering of nodes
    for (size_t i=0; i<size(FreeSurface_Nodes); i++) FreeSurface_Nodes[i]->ID = i;

    if (Panels_Outward){
        for (SP_Geo G : FreeSurface_Elements) G->Reorder_Nodes();
    }

    // for (int i=0; i<size(FreeSurface_Nodes); i++) std::cout <<     FreeSurface_Nodes[i]->ID << std::endl;
}

}
