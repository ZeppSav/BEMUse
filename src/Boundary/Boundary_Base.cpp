//---------------------------------------------------------------
//------------------ Boundary Functions--------------------------
//---------------------------------------------------------------

#include "Boundary_Base.h"

namespace BEMUse
{

//--- Parameter template functions
template <> int     Parameter::Get_Param<int>()    {return std::get<1>(*this);}
template <> Real    Parameter::Get_Param<Real>()   {return std::get<2>(*this);}
template <> bool    Parameter::Get_Param<bool>()   {return std::get<3>(*this);}

//--- Geo setup

void Boundary::Setup()
{
    // This is the generic function to prepare the geometry

    //--- Calculate geometric & kinematic quantities
    Generate_Nodes();
    Generate_Elements();
    Generate_Aux_Nodes();
    Generate_Aux_Elements();
    Generate_FreeSurface_Nodes();
    Generate_FreeSurface_Elements();
    Unify_Geometries();

    //--- Calculate geometric & kinematic quantities
    for (SP_Geo G : Elements)       G->Set_Centroid();
    for (SP_Geo G : Aux_Elements)   G->Set_Centroid();
    for (SP_Geo G : FreeSurface_Elements)   G->Set_Centroid();

    Calculate_GridParams();
    Calculate_Volume();
    Calculate_SurfaceArea();
    Calculate_COB();
    Set_Kin_Params();

    std::cout   << "Geometry generated. "
                << Elements.size() << " elements, "
                << Aux_Elements.size() << " auxiliary (free surface) elements. "
                << FreeSurface_Elements.size() << " external (free surface) elements. \n";
    std::cout   << "Surface area = " << Surface_Area
                << " m^2. Volume = " << Volume << " m^3"    << std::endl;
}

//--- Geometry functions

void Boundary::Generate_Symmetry_Nodes(PLANE P)
{
    // The geometry having been defined in the initial steps, nodes are reflected about a plane.

    for (SP_Node N : Nodes){
        Vector3 PR = N->Position_Global();
        if (P==XPlane)  PR(0) *= -1.0;
        if (P==YPlane)  PR(1) *= -1.0;
        if (P==ZPlane)  PR(2) *= -1.0;
        Symm_Nodes.push_back(std::make_shared<Node>(Global_CS,PR));
    }
}

void Boundary::Generate_Symmetry_Elements()
{
    // The geometry having been defined in the initial steps, elements are reflected about a plane.
    // The only application case as this stage is hydrodynamics, so the elements shall be reflected about the plane z=0
    // Nb. Node order must be reversed to ensure normal is facing in opposite direction.

    for (SP_Geo G : Elements){
        if (G->Get_N()==3){
            SP_Node RNode1 = Symm_Nodes[G->Get_Node(2)->ID];
            SP_Node RNode2 = Symm_Nodes[G->Get_Node(1)->ID];
            SP_Node RNode3 = Symm_Nodes[G->Get_Node(0)->ID];
            Symm_Elements.push_back(std::make_shared<Tri_Element>(RNode1,RNode2,RNode3));
        }
        if (G->Get_N()==4){
            SP_Node RNode1 = Symm_Nodes[G->Get_Node(3)->ID];
            SP_Node RNode2 = Symm_Nodes[G->Get_Node(2)->ID];
            SP_Node RNode3 = Symm_Nodes[G->Get_Node(1)->ID];
            SP_Node RNode4 = Symm_Nodes[G->Get_Node(0)->ID];
            Symm_Elements.push_back(std::make_shared<Quad_Element>(RNode1,RNode2,RNode3,RNode4));
        }
    }
}

void Boundary::Generate_Symmetry_Aux_Nodes(PLANE P)
{
    // The geometry having been defined in the initial steps, nodes are reflected about a plane.

    for (SP_Node N : Aux_Nodes){
        Vector3 PR = N->Position_Global();
        if (P==XPlane)  PR(0) *= -1.0;
        if (P==YPlane)  PR(1) *= -1.0;
        if (P==ZPlane)  PR(2) *= -1.0;
        Symm_Aux_Nodes.push_back(std::make_shared<Node>(Global_CS,PR));
    }
}

void Boundary::Generate_Symmetry_Aux_Elements()
{
    // The geometry having been defined in the initial steps, elements are reflected about a plane.
    // The only application case as this stage is hydrodynamics, so the elements shall be reflected about the plane z=0
    // Nb. Node order must be reversed to ensure normal is facing in opposite direction.

    for (SP_Geo G : Aux_Elements){
        if (G->Get_N()==3){
            SP_Node RNode1 = Symm_Aux_Nodes[G->Get_Node(2)->ID];
            SP_Node RNode2 = Symm_Aux_Nodes[G->Get_Node(1)->ID];
            SP_Node RNode3 = Symm_Aux_Nodes[G->Get_Node(0)->ID];
            Symm_Aux_Elements.push_back(std::make_shared<Tri_Element>(RNode1,RNode2,RNode3));
        }
        if (G->Get_N()==4){
            SP_Node RNode1 = Symm_Aux_Nodes[G->Get_Node(3)->ID];
            SP_Node RNode2 = Symm_Aux_Nodes[G->Get_Node(2)->ID];
            SP_Node RNode3 = Symm_Aux_Nodes[G->Get_Node(1)->ID];
            SP_Node RNode4 = Symm_Aux_Nodes[G->Get_Node(0)->ID];
            Symm_Aux_Elements.push_back(std::make_shared<Quad_Element>(RNode1,RNode2,RNode3,RNode4));
        }
    }
}

//--- Geometry calculations

void Boundary::Calculate_GridParams()
{
    // This calculates various important geometric parameters of the body

    // What are the largest spatial dimensions? This allows for a global size parameter

    StdVector X,Y,Z;
    for (SP_Node N : Nodes){
        Vector3 P = N->Position_Global();
        X.push_back(P(0));
        Y.push_back(P(1));
        Z.push_back(P(2));
    }

    Real Xmin = *std::min_element(X.begin(), X.end());
    Real Xmax = *std::max_element(X.begin(), X.end());
    Real Ymin = *std::min_element(Y.begin(), Y.end());
    Real Ymax = *std::max_element(Y.begin(), Y.end());
    Real Zmin = *std::min_element(Z.begin(), Z.end());
    Real Zmax = *std::max_element(Z.begin(), Z.end());

    // Set grid maximum dimension
    if ((Xmax-Xmin)>Max_Dim) Max_Dim = (Xmax-Xmin);
    if ((Ymax-Ymin)>Max_Dim) Max_Dim = (Ymax-Ymin);
    if ((Zmax-Zmin)>Max_Dim) Max_Dim = (Zmax-Zmin);

    // Now repeat the same for auxiliary elements
    if (Aux_Nodes.empty())  return;

    StdVector AX,AY,AZ;
    for (SP_Node N : Aux_Nodes){
        Vector3 P = N->Position_Global();
        AX.push_back(P(0));
        AY.push_back(P(1));
        AZ.push_back(P(2));
    }

    Real AXmin = *std::min_element(AX.begin(), AX.end());
    Real AXmax = *std::max_element(AX.begin(), AX.end());
    Real AYmin = *std::min_element(AY.begin(), AY.end());
    Real AYmax = *std::max_element(AY.begin(), AY.end());
    Real AZmin = *std::min_element(AZ.begin(), AZ.end());
    Real AZmax = *std::max_element(AZ.begin(), AZ.end());

    // Set grid maximum dimension
    if ((AXmax-AXmin)>Max_Aux_Dim) Max_Aux_Dim = (Xmax-Xmin);
    if ((AYmax-AYmin)>Max_Aux_Dim) Max_Aux_Dim = (Ymax-Ymin);
    if ((AZmax-AZmin)>Max_Aux_Dim) Max_Aux_Dim = (Zmax-Zmin);

//    std::cout   << "Maximum volume dimension = " << Max_Dim << " Maximum auxiliary dimension = " << Max_Aux_Dim << ".\n";
}

void Boundary::Calculate_Volume()
{
    // What is the volume enclosed by the boundary? This is a simple surface integral

    // Calculate volume for each normal direction
    Real Vol_x = 0, Vol_y = 0, Vol_z = 0;

    for (SP_Geo G : Elements)
    {
        Vector3 Normal = G->Centroid->Z_Axis_Global();
        Vector3 Pos = G->Centroid->Position_Global();
        Real A = G->Area;

        // Volume
        Vol_x += -Pos(0)*Normal(0)*A;
        Vol_y += -Pos(1)*Normal(1)*A;
        Vol_z += -Pos(2)*Normal(2)*A;
    }

    //--- Volume
    Volume = (Vol_x+Vol_y+Vol_z)/3.0;
}

void Boundary::Calculate_COB()
{
    // What is the centre of buoyancy of the body?

    for (SP_Geo G : Elements)
    {
        Vector3 Normal = G->Centroid->Z_Axis_Global();
        Vector3 Pos = G->Centroid->Position_Global();
        Real A = G->Area;

        // Centre of buoyancy
        Centre_Buoyancy(0) -= Normal(0)*Pos(0)*Pos(0)*A;
        Centre_Buoyancy(1) -= Normal(1)*Pos(1)*Pos(1)*A;
        Centre_Buoyancy(2) -= Normal(2)*Pos(2)*Pos(2)*A;
    }

    Centre_Buoyancy *= 0.5/Volume;

    // Clean up numerical zeros
    if (fabs(Centre_Buoyancy(0))<Max_Dim*1e-6)  Centre_Buoyancy(0) = 0;
    if (fabs(Centre_Buoyancy(1))<Max_Dim*1e-6)  Centre_Buoyancy(1) = 0;
    if (fabs(Centre_Buoyancy(2))<Max_Dim*1e-6)  Centre_Buoyancy(2) = 0;

    // For now I will assume that the centre of gravity has either not been set or has been ignored
}

void Boundary::Calculate_SurfaceArea()
{
    // What is the surface area of the boundary?

    // Sum area
    for (SP_Geo G : Elements) Surface_Area += G->Area;
}

//--- Import/Export

void Boundary::Export_Geometry_GDF(std::string &FilePath)
{
    // This exports the geometry in .gdf (WAMIT) format

//    std::string Path = FilePath + ".bem";
    std::ofstream file;

    file.open(FilePath, std::ofstream::out | std::ofstream::trunc); // Clear!
    file.close();

    // Write initial inputs
    file.open(FilePath);
    if (file.is_open())
    {
        file << "Geometry file .gdf type for use in WAMIT flow solver. The file was generated in BEMUse. Viva la BEMUse!\n";
        file << "1 9.80665 	ULEN GRAV\n";
        file << "0  0 	ISX  ISY\n";      // Assume no symmetries
        file << Elements.size() << std::endl;
        for (SP_Geo G : Elements)
        {
            StateVector N;
            G->Get_Nodal_Positions(N);
            if (G->Get_N()==3){
                file << N[0](0) << "  " << N[0](1) << "  " << N[0](2) << std::endl;
                file << N[1](0) << "  " << N[1](1) << "  " << N[1](2) << std::endl;
                file << N[2](0) << "  " << N[2](1) << "  " << N[2](2) << std::endl;
                file << N[2](0) << "  " << N[2](1) << "  " << N[2](2) << std::endl;
            }
            if (G->Get_N()==4){
                file << N[0](0) << "  " << N[0](1) << "  " << N[0](2) << std::endl;
                file << N[1](0) << "  " << N[1](1) << "  " << N[1](2) << std::endl;
                file << N[2](0) << "  " << N[2](1) << "  " << N[2](2) << std::endl;
                file << N[3](0) << "  " << N[3](1) << "  " << N[3](2) << std::endl;
            }
        }
    }
    file.close();
}

void Boundary::Export_Geometry_MAR(std::string &FilePath)
{
    // This exports the geometry in .mar (NEMOH) format

    std::ofstream file;

    file.open(FilePath, std::ofstream::out | std::ofstream::trunc); // Clear!
    file.close();

    // Write initial inputs
    file.open(FilePath);
    if (file.is_open())
    {
        file << "                    2          0\n";       // Symmetry!!!
        for (SP_Node N : Nodes){
            Vector3 P = N->Position_Global();
            file  << N->ID+1 << "  " << P(0) << "  " << P(1) << "  " << P(2) << std::endl;
        }
        file << "             0          0.00          0.00          0.00\n";
        for (SP_Geo G : Elements){
            if (G->Get_N()==3) file     << G->Get_Node(0)->ID+1 << "  " \
                                        << G->Get_Node(1)->ID+1 << "  " \
                                        << G->Get_Node(2)->ID+1 << "  " \
                                        << G->Get_Node(2)->ID+1 << std::endl;
            if (G->Get_N()==4) file     << G->Get_Node(0)->ID+1 << "  " \
                                        << G->Get_Node(1)->ID+1 << "  " \
                                        << G->Get_Node(2)->ID+1 << "  " \
                                        << G->Get_Node(3)->ID+1 << std::endl;
        }
        file << "          0          0          0          0\n";
    }
    file.close();
}

//--- Destructor

Boundary::~Boundary()
{
    // Ad advantage of using shared pointer for all objects is that we simply need to remove the reference to the object
    // in order to delete it.

    // Clear elements
    Elements.clear();
    Aux_Elements.clear();
    Symm_Elements.clear();
    FreeSurface_Elements.clear();

    // Clear nodes
    Nodes.clear();
    Aux_Nodes.clear();
    Symm_Nodes.clear();
    FreeSurface_Nodes.clear();
}

}
