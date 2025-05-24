//-----------------------------------------------------------------------------
//-------------------------PNL Geo routines------------------------------------
//-----------------------------------------------------------------------------

#include "PNL_Geo.h"

namespace BEMUse
{

//--- Import/Export

inline bool contains(const std::string& source, const std::string& target) {
    return source.find(target) != std::string::npos;
}

bool PNL_Geometry::Read_Input_File(std::string &FilePath)
{
    // This function opens the File with filename in the compile directory and reads
    // through the file storing node positions and panel connectivity

    PNL_Nodes_Raw.clear();
    PNL_Connectivity.clear();

    std::ifstream file(FilePath);
    std::string line;
    int Pos = 0;
    if (file.is_open())
    {
        while ( std::getline (file,line) )
        {
            // Split line up
            if (line.empty()) continue;
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments

            // Check & mark positions
            if (contains(line,"Number of Panels")){
                // The next line contains the number of panels, nodes & X & Y Symmetry
                Pos = 1;
                continue;
            }
            if (contains(line,"Start Definition of Node Coordinates")){
                // The next line contains the number of panels, nodes & X & Y Symmetry
                Pos = 2;
                continue;
            }
            if (contains(line,"End Definition of Node Coordinates")){
                // The next line contains the number of panels, nodes & X & Y Symmetry
                Pos = 3;
                continue;
            }
            if (contains(line,"Start Definition of Node Relations")){
                // The next line contains the number of panels, nodes & X & Y Symmetry
                Pos = 4;
                continue;
            }
            if (contains(line,"End Definition of Node Relations")){
                // The next line contains the number of panels, nodes & X & Y Symmetry
                Pos = 5;
                continue;
            }

            // Line contains the number of panels, nodes & X & Y Symmetry
            if (Pos == 1){
                NPanels = std::stoi(Fields[0]);
                NNodes = std::stoi(Fields[1]);
                if (std::stoi(Fields[2])) XSymmetry = true;
                if (std::stoi(Fields[3])) YSymmetry = true;
            }

            // Node coordinate
            if (Pos == 2){
                Vector3 Pos(std::stod(Fields[1]),std::stod(Fields[2]),std::stod(Fields[3]));
                PNL_Nodes_Raw.push_back(Pos);
            }

            // Panel connectivity
            if (Pos == 4){
                std::vector<int> NDS;
                NDS.push_back(std::stoi(Fields[2]));
                NDS.push_back(std::stoi(Fields[3]));
                NDS.push_back(std::stoi(Fields[4]));
                if (std::stoi(Fields[1])==4)    NDS.push_back(std::stoi(Fields[5]));
                PNL_Connectivity.push_back(NDS);
            }
        }
        file.close();
        return true;
    }
    else{
        std::cout << FilePath << " . A suitable .pnl file was not found at this position.\n";
        return false;
    }
}

//--- Geometry generation

void PNL_Geometry::Generate_Nodes()
{
    // This stores the nodes of the pnl import

    Inertial_CS = new CoordSys();

    Real Scale = 1.0;
    for (int i=0; i<PNL_Nodes_Raw.size(); i++){
        Vector3 PS = Scale*PNL_Nodes_Raw[i];
        Nodes.push_back(std::make_shared<Node>(Inertial_CS,PS));
    }
    for (size_t i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

    PNL_Nodes_Raw.clear();        // Clear list
}

void PNL_Geometry::Generate_Elements()
{
    // This generates panel elements from the pnl nodes

    for (int i=0; i<PNL_Connectivity.size(); i++){

        if (PNL_Connectivity[i].size()==3)   Elements.push_back(std::make_shared<Tri_Element>(  Nodes[PNL_Connectivity[i][0]-1],
                                                                                                Nodes[PNL_Connectivity[i][1]-1],
                                                                                                Nodes[PNL_Connectivity[i][2]-1]));
        if (PNL_Connectivity[i].size()==4)   Elements.push_back(std::make_shared<Quad_Element>( Nodes[PNL_Connectivity[i][0]-1],
                                                                                                Nodes[PNL_Connectivity[i][1]-1],
                                                                                                Nodes[PNL_Connectivity[i][2]-1],
                                                                                                Nodes[PNL_Connectivity[i][3]-1]));
    }

    PNL_Connectivity.clear();            // Clear list to avoid problems with destructor
}

}
