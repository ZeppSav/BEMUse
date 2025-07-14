//-----------------------------------------------------------------------------
//-------------------------Wing Class------------------------------------------
//-----------------------------------------------------------------------------

#include "Wing.h"

namespace BEMUse
{

// Parameter specification

void Wing::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Span"))             Span = P.Get_Param<Real>();
        if (P.myNameis("BaselineChord"))    C = P.Get_Param<Real>();
        if (P.myNameis("Cosine_Disc"))      Cosine = P.Get_Param<bool>();
        if (P.myNameis("NPanels_Span"))     NS = P.Get_Param<int>();
        if (P.myNameis("NPanels_Chord"))    NC = P.Get_Param<int>();
    }
}

//--- Geometry functions

void Wing::SetSpan()
{
    // This sets the spanwise positions of the profiles
    if (Cosine) for (int i=0; i<NS+1; i++) SpanPos.push_back(-Span+2.0*Span*PiCosFac(i,NS+1));
    else        for (int i=0; i<NS+1; i++) SpanPos.push_back(-Span+2.0*Span*i/NS);
}

void Wing::SetChord()
{
    // This sets the chordwise distribution as a function of span position

//    for (int i=0; i<NS+1; i++) Chord.push_back(1.0);  // Rectangular
    for (int i=0; i<NS; i++){
        Real xnorm = SpanPos[i]/(Span);
        Chord.push_back(sqrt(C*(1-xnorm*xnorm)));         // Elliptical
    }

}

void Wing::SetAlpha0()
{
    // This set the geometric angle of attack of the airfoil

    for (int i=0; i<NS+1; i++) Alpha0.push_back(0.0);  // Zero
}

//--- Geometry generation

void Wing::Generate_Nodes()
{
    // This generates the nodes of the airfoil surface. We however first need to generate the sectional data.

    SetSpan();
    SetChord();
    SetAlpha0();

    if (!Global_CS)    Global_CS = new CoordSys();
//    Quat Q = Quat::Identity();
    Quat Q = Quat::FromTwoVectors(UnitX,Vector3(cos(-15*D2R),0.0,sin(-15*D2R)));
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS,Q,Vector3::Zero());

    CReal C_KT(-0.4,0.1);     // Centre for a karman trefftz profile
    Real TE_A = 15*D2R;
    Real Shift = 0.25;      // Factor to shift airfoil in x direction

    for (int i=0; i<NS+1; i++)
    {
        // Generate for each spanwise section the corresponding coordinate system and airfoil coordinates.

        // Generate profile
        StdVector X,Z;
        Karman_Trefftz(C_KT,TE_A,NC,X,Z);        // Karman Trefftz profile

        // Scale profile
        Normalize(NC+1,X,Z);
        for (int j=0; j<NC+1; j++) X[j] -= Shift;
        for (int j=0; j<NC+1; j++) X[j] *= Chord[i];
        for (int j=0; j<NC+1; j++) Z[j] *= Chord[i];

        // Generate nodes
//        for (int j=0; j<NC+1; j++) Nodes.push_back(std::make_shared<Node>(CSLoc,Coords[j]));
        for (int j=0; j<NC+1; j++)  Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(X[j],SpanPos[i],Z[j])));
    }

    // Set node ID
    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

}

void Wing::Generate_Elements()
{
    // This generates the elements of the airfoil surface.

    for (int i=0; i<NS; i++){
        for (int j=0; j<NC; j++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(j,i)],
                                                                                        Nodes[Node_ID(j,i+1)],
                                                                                        Nodes[Node_ID(j+1,i+1)],
                                                                                        Nodes[Node_ID(j+1,i)]));
    }


}

}
