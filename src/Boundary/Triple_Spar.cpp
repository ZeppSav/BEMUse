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
void Triple_Spar::Set_Discretisation(std::vector<int> &D)
{
    // Specifies discretisation of the different elements
    Leg1->Set_Discretisation(D);
    Leg2->Set_Discretisation(D);
    Leg3->Set_Discretisation(D);
}

void Triple_Spar::Set_Auxiliary_Discretisation(std::vector<int> &D)
{
    // The sections are all volumes of revolution, so this simply specifies for all bodies the internal radial discretisation
    if (Central_Cylinder) Central_Cylinder->Set_Auxiliary_Discretisation(D);
    Leg1->Set_Auxiliary_Discretisation(D);
    Leg2->Set_Auxiliary_Discretisation(D);
    Leg3->Set_Auxiliary_Discretisation(D);
}

void Triple_Spar::Set_External_Discretisation(std::vector<int> &D)
{
    // The individual external free surface discretisations are not used.
    // This is generated differently for the OC4 platform (must be implemented)

    NAE = D[0]; NRE = D[1];
}

void Triple_Spar::Set_Dimensions(std::vector<Real> &D)
{
    // Dimensions

    Leg1->Set_Dimensions(D);
    Leg2->Set_Dimensions(D);
    Leg3->Set_Dimensions(D);

    // Set leg centres
    R_Leg = D[2];
    Real Theta1 = -0.5*PI, Theta2 = -0.5*PI+2.0*PI/3, Theta3 = -0.5*PI+4.0*PI/3;
    C1 = Vector3(R_Leg*sin(Theta1),R_Leg*cos(Theta1),0.0);
    C2 = Vector3(R_Leg*sin(Theta2),R_Leg*cos(Theta2),0.0);
    C3 = Vector3(R_Leg*sin(Theta3),R_Leg*cos(Theta3),0.0);

    // This flag shifts the position of the floater such that the origin is directly above the first leg.
    // Some turbine floater designs are constructed like this to provide stiffness in pitch direction (among other reasons).
    if (ShiftRearwards){
        C1 -= Vector3(R_Leg*sin(Theta1),0.0,0.0);
        C2 -= Vector3(R_Leg*sin(Theta1),0.0,0.0);
        C3 -= Vector3(R_Leg*sin(Theta1),0.0,0.0);
    }
}

void Triple_Spar::Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0];}

void Triple_Spar::Set_Flags(std::vector<bool> &D)
{
    ShiftRearwards = D[1];
    if (Central_Cylinder) Central_Cylinder->Set_Flags(D);
    Leg1->Set_Flags(D);
    Leg2->Set_Flags(D);
    Leg3->Set_Flags(D);
}

//--- Geometry functions
void Triple_Spar::Generate_Nodes()
{
    // Generate the nodes within each geometry

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);
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
void SemiSub_OC4::Set_Discretisation(std::vector<int> &D)
{
    // Specifies discretisation of the different elements

//    NA = D[0];

    std::vector<int> CentreDisc;
    CentreDisc.push_back(D[0]);
    CentreDisc.push_back(D[1]);
    CentreDisc.push_back(D[4]);
    Central_Cylinder->Set_Discretisation(CentreDisc);

    std::vector<int> PlanetDisc;
    PlanetDisc.push_back(D[0]);
    PlanetDisc.push_back(D[1]);
    PlanetDisc.push_back(D[2]);
    PlanetDisc.push_back(D[3]);
    PlanetDisc.push_back(D[4]);

    Leg1->Set_Discretisation(PlanetDisc);
    Leg2->Set_Discretisation(PlanetDisc);
    Leg3->Set_Discretisation(PlanetDisc);
}

void SemiSub_OC4::Set_Dimensions(std::vector<Real> &D)
{
    // Dimensions
    ScaleFac = D[0];

    std::vector<Real> CentreDims, PlanetDims;
    CentreDims.push_back(3.25*ScaleFac);
    CentreDims.push_back(20*ScaleFac);

    PlanetDims.push_back(12*ScaleFac);
    PlanetDims.push_back(6*ScaleFac);
    PlanetDims.push_back(20*ScaleFac);
    PlanetDims.push_back(6*ScaleFac);

    Central_Cylinder->Set_Dimensions(CentreDims);
    Leg1->Set_Dimensions(PlanetDims);
    Leg2->Set_Dimensions(PlanetDims);
    Leg3->Set_Dimensions(PlanetDims);

    // Set leg centres
    R_Leg = 28.87*ScaleFac;
    Real Theta1 = -0.5*PI, Theta2 = -0.5*PI+2.0*PI/3, Theta3 = -0.5*PI+4.0*PI/3;
    C1 = Vector3(R_Leg*sin(Theta1),R_Leg*cos(Theta1),0.0);
    C2 = Vector3(R_Leg*sin(Theta2),R_Leg*cos(Theta2),0.0);
    C3 = Vector3(R_Leg*sin(Theta3),R_Leg*cos(Theta3),0.0);
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
