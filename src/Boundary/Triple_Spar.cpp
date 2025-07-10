//-----------------------------------------------------------------------------
//-------------------------Triple Spar routines--------------------------------
//-----------------------------------------------------------------------------

#include "Triple_Spar.h"

namespace BEMUse
{

//--- Constructor

Triple_Spar::Triple_Spar()
{
    // Generate geometry objects
    Leg1 = new Half_Cylinder();
    Leg2 = new Half_Cylinder();
    Leg3 = new Half_Cylinder();
}

//--- Geometry specification

void Triple_Spar::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    if (Central_Cylinder) Central_Cylinder->Set_Parameters(Params);
    Leg1->Set_Parameters(Params);
    Leg2->Set_Parameters(Params);
    Leg3->Set_Parameters(Params);

    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Radius"))           R_Loc = P.Get_Param<Real>();
        if (P.myNameis("Spar_Radius"))      R_Leg = P.Get_Param<Real>();
        if (P.myNameis("Shift_Over_Leg"))   ShiftRearwards = P.Get_Param<bool>();
    }
}

//--- Geometry functions

void Triple_Spar::Generate_Nodes()
{
    // Generate the nodes within each geometry

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    // Set leg centres
    Real Theta1 = -0.5*PI, Theta2 = -0.5*PI+2.0*PI/3, Theta3 = -0.5*PI+4.0*PI/3;
    C1 = Vector3(R_Leg*sin(Theta1),R_Leg*cos(Theta1),0.0);
    C2 = Vector3(R_Leg*sin(Theta2),R_Leg*cos(Theta2),0.0);
    C3 = Vector3(R_Leg*sin(Theta3),R_Leg*cos(Theta3),0.0);

    if (ShiftRearwards){
        C1 -= Vector3(R_Leg*sin(Theta1),0.0,0.0);
        C2 -= Vector3(R_Leg*sin(Theta1),0.0,0.0);
        C3 -= Vector3(R_Leg*sin(Theta1),0.0,0.0);
    }

    CSLeg1 = new CoordSys(Global_CS,Quat::Identity(),C1);
    CSLeg2 = new CoordSys(Global_CS,Quat::Identity(),C2);
    CSLeg3 = new CoordSys(Global_CS,Quat::Identity(),C3);

    if (Central_Cylinder) Central_Cylinder->Set_CoordinateSystems(Global_CS, Inertial_CS);
    Leg1->Set_CoordinateSystems(Global_CS, CSLeg1);
    Leg2->Set_CoordinateSystems(Global_CS, CSLeg2);
    Leg3->Set_CoordinateSystems(Global_CS, CSLeg3);

    // Generate Nodes
    if (Central_Cylinder) Central_Cylinder->Generate_Nodes();
    Leg1->Generate_Nodes();
    Leg2->Generate_Nodes();
    Leg3->Generate_Nodes();
}

void Triple_Spar::Generate_Aux_Nodes()
{
    // Generate auxiliary Nodes
    if (Central_Cylinder) Central_Cylinder->Generate_Aux_Nodes();
    Leg1->Generate_Aux_Nodes();
    Leg2->Generate_Aux_Nodes();
    Leg3->Generate_Aux_Nodes();
}

void Triple_Spar::Generate_Elements()
{
    if (Central_Cylinder) Central_Cylinder->Generate_Elements();
    Leg1->Generate_Elements();
    Leg2->Generate_Elements();
    Leg3->Generate_Elements();
}

void Triple_Spar::Generate_Aux_Elements()
{
    if (Central_Cylinder) Central_Cylinder->Generate_Aux_Elements();
    Leg1->Generate_Aux_Elements();
    Leg2->Generate_Aux_Elements();
    Leg3->Generate_Aux_Elements();
}

void Triple_Spar::Generate_FreeSurface_Nodes()
{
    // Generate auxiliary Nodes
    if (Central_Cylinder) Central_Cylinder->Generate_FreeSurface_Nodes();
    Leg1->Generate_FreeSurface_Nodes();
    Leg2->Generate_FreeSurface_Nodes();
    Leg3->Generate_FreeSurface_Nodes();
}

void Triple_Spar::Generate_FreeSurface_Elements()
{
    if (Central_Cylinder) Central_Cylinder->Generate_FreeSurface_Elements();
    Leg1->Generate_FreeSurface_Elements();
    Leg2->Generate_FreeSurface_Elements();
    Leg3->Generate_FreeSurface_Elements();
}

//--- Unify multiple geos

void Triple_Spar::Unify_Geometries()
{
    // Nodes
    if (Central_Cylinder) Central_Cylinder->Get_Nodes(Nodes);
    Leg1->Get_Nodes(Nodes);
    Leg2->Get_Nodes(Nodes);
    Leg3->Get_Nodes(Nodes);
    for (int i=0; i<Nodes.size(); i++)      Nodes[i]->ID = i;       // Reset for panel deflection

    // Aux nodes
    if (Central_Cylinder) Central_Cylinder->Get_Aux_Nodes(Aux_Nodes);
    Leg1->Get_Aux_Nodes(Aux_Nodes);
    Leg2->Get_Aux_Nodes(Aux_Nodes);
    Leg3->Get_Aux_Nodes(Aux_Nodes);
    for (int i=0; i<Aux_Nodes.size(); i++)  Aux_Nodes[i]->ID = i;

    // Ext nodes
    if (Central_Cylinder) Central_Cylinder->Get_FreeSurface_Nodes(FreeSurface_Nodes);
    Leg1->Get_FreeSurface_Nodes(FreeSurface_Nodes);
    Leg2->Get_FreeSurface_Nodes(FreeSurface_Nodes);
    Leg3->Get_FreeSurface_Nodes(FreeSurface_Nodes);
    for (int i=0; i<FreeSurface_Nodes.size(); i++)  FreeSurface_Nodes[i]->ID = i;

    // Elements
    if (Central_Cylinder) Central_Cylinder->Get_Elements(Elements);
    Leg1->Get_Elements(Elements);
    Leg2->Get_Elements(Elements);
    Leg3->Get_Elements(Elements);

    // Aux Elements
    if (Central_Cylinder) Central_Cylinder->Get_Aux_Elements(Aux_Elements);
    Leg1->Get_Aux_Elements(Aux_Elements);
    Leg2->Get_Aux_Elements(Aux_Elements);
    Leg3->Get_Aux_Elements(Aux_Elements);

    // Aux Elements
    if (Central_Cylinder) Central_Cylinder->Get_FreeSurface_Elements(FreeSurface_Elements);
    Leg1->Get_FreeSurface_Elements(FreeSurface_Elements);
    Leg2->Get_FreeSurface_Elements(FreeSurface_Elements);
    Leg3->Get_FreeSurface_Elements(FreeSurface_Elements);

    // We can delete these objects now as they are no longer necessary
    if (Central_Cylinder) delete Central_Cylinder;
    delete Leg1;
    delete Leg2;
    delete Leg3;
}

//--- Constructor

SemiSub_OC4::SemiSub_OC4()
{
    // Generate geometry objects
    Central_Cylinder = new Half_Cylinder();
    Leg1 = new Spar_Leg();
    Leg2 = new Spar_Leg();
    Leg3 = new Spar_Leg();
}

//--- Geometry specification

void SemiSub_OC4::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    // We import the scalefac, everything else is known.

    for (Parameter P : Parameters)
    {
        if (P.myNameis("Scaling_Factor"))   ScaleFac = P.Get_Param<Real>();
    }

    // Parameters central column
    std::vector<Parameter> CentralParams;
    StdAppend(CentralParams, Params);
    CentralParams.push_back(Parameter("Radius", ScaleFac*Real(3.25)));
    CentralParams.push_back(Parameter("Draft", ScaleFac*Real(20.0)));
    Central_Cylinder->Set_Parameters(CentralParams);

    // Set spar leg parameters
    std::vector<Parameter> SparLegParams;
    StdAppend(SparLegParams, Params);
    SparLegParams.push_back(Parameter("Leg_Lower_Radius", ScaleFac*Real(12.0)));
    SparLegParams.push_back(Parameter("Leg_Upper_Radius", ScaleFac*Real(6.0)));
    SparLegParams.push_back(Parameter("Draft", ScaleFac*Real(20.0)));
    SparLegParams.push_back(Parameter("Heave_Plate_Height", ScaleFac*Real(6.0)));
    Leg1->Set_Parameters(SparLegParams);
    Leg2->Set_Parameters(SparLegParams);
    Leg3->Set_Parameters(SparLegParams);

    // Set leg centres
    R_Leg = 28.87*ScaleFac;
}

//--- Geometry calculations (Hard coded)

void SemiSub_OC4::Set_Kin_Params()
{
    // Hard-code here the kinematic parameters (Note: These are configured for density of 1025kg/m^3.)
    Centre_Gravity = Vector3(0,0,-13.46);
    Mass = 1.34730e+07;

    Mass_Mat(0,0) = Mass;
    Mass_Mat(1,1) = Mass;
    Mass_Mat(2,2) = Mass;
    Mass_Mat(3,3) = 6.82700e+09;
    Mass_Mat(4,4) = 6.82700e+09;
    Mass_Mat(5,5) = 1.22600e+10;

    Stiffness_Mat(2,2) = 3.836e+06 ;
    Stiffness_Mat(3,3) = -3.776e+08;
    Stiffness_Mat(4,4) = -3.776e+08;

    // Kinemetic scaling
    Real SF3 = ScaleFac*ScaleFac*ScaleFac;      // Volume scaling
    Real SF5 = SF3*ScaleFac*ScaleFac;           // Moment of inertia scaling

    Centre_Gravity(2) *= ScaleFac;
    Mass *= SF3;
    Mass_Mat.block(0,0,3,3) *= SF3;
    Mass_Mat.block(3,3,3,3) *= SF5;
    Stiffness_Mat.block(0,0,3,3) *= SF3;
    Stiffness_Mat.block(3,3,3,3) *= SF5;
}

}
