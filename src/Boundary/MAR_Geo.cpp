//-----------------------------------------------------------------------------
//-------------------------MAR Geo routines------------------------------------
//-----------------------------------------------------------------------------

#include "MAR_Geo.h"

namespace BEMUse
{

//--- Import/Export

bool MAR_Geometry::Read_Input_File(std::string &FilePath)
{
    // This function opens the File with filename in the compile directory and reads
    // through the file storing node positions and panel connectivity

    MAR_Nodes_Raw.clear();
    MAR_Connectivity.clear();

    std::ifstream file(FilePath);
    std::string line;
    if (file.is_open())
    {
        bool Connectivity = false;
        while ( std::getline (file,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields.size()==2)    continue;  // Read symmetrix info.
            if (Fields.size()==4){
                int ID = std::stoi(Fields[0]);
                if (ID==0){
                    Connectivity=true;
                    continue;
                }

                if (Connectivity) {
                    std::vector<int> NDS;
                    NDS.push_back(std::stoi(Fields[0])-1);     // Index 1!
                    NDS.push_back(std::stoi(Fields[1])-1);
                    NDS.push_back(std::stoi(Fields[2])-1);
                    if (Fields[2]!=Fields[3]) NDS.push_back(std::stoi(Fields[3])-1);
                    MAR_Connectivity.push_back(NDS);
                }
                else {
                    Vector3 Pos(std::stod(Fields[1]),std::stod(Fields[2]),std::stod(Fields[3]));
                    MAR_Nodes_Raw.push_back(Pos);
                }
            }

        }
        file.close();
        return true;
    }
    else{
        std::cout << FilePath << " . A suitable .mar file was not found at this position.\n";
        return false;
    }

}

//--- Geometry generation

void MAR_Geometry::Generate_Nodes()
{
    // This temporarily stores the nodes of the MAR import

    Inertial_CS = new CoordSys();

    Real Scale = 1.0;
    for (int i=0; i<MAR_Nodes_Raw.size(); i++){
        Vector3 PS = Scale*MAR_Nodes_Raw[i];
        Nodes.push_back(std::make_shared<Node>(Inertial_CS,PS));
    }
    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

    MAR_Nodes_Raw.clear();        // Clear list
}

void MAR_Geometry::Generate_Elements()
{
    // This generates panel elements from the MAR nodes

    for (int i=0; i<MAR_Connectivity.size(); i++){

        if (MAR_Connectivity[i].size()==3)   Elements.push_back(std::make_shared<Tri_Element>(  Nodes[MAR_Connectivity[i][0]],
                                                                                                Nodes[MAR_Connectivity[i][1]],
                                                                                                Nodes[MAR_Connectivity[i][2]]));
        if (MAR_Connectivity[i].size()==4)   Elements.push_back(std::make_shared<Quad_Element>( Nodes[MAR_Connectivity[i][0]],
                                                                                                Nodes[MAR_Connectivity[i][1]],
                                                                                                Nodes[MAR_Connectivity[i][2]],
                                                                                                Nodes[MAR_Connectivity[i][3]]));
    }

    MAR_Connectivity.clear();            // Clear list to avoid problems with destructor
}

}
