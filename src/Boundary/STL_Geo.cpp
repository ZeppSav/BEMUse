//-----------------------------------------------------------------------------
//-------------------------STL Geo routines------------------------------------
//-----------------------------------------------------------------------------

#include "STL_Geo.h"

namespace BEMUse
{

//--- Import/Export

bool STL_Geometry::Read_Input_File(std::string &FilePath)
{
    // This function opens the File with filename in the compile directory and reads
    // through the file storing node positions and panel connectivity

    STL_Els_Raw.clear();
    STL_Els.clear();

    std::ifstream file(FilePath);
    std::string line;
    if (file.is_open())
    {

        int count;
        int countglob = 0;
        std::vector<Vector3> El;

        while ( std::getline (file,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields[0]=="facet")
            {
                count=0;
        //                El.push_back(Vector3(std::stod(Fields[2]),std::stod(Fields[3]),std::stod(Fields[4])));      // Normal
            }
            if (count==2)   El.push_back(Vector3(std::stod(Fields[1]),std::stod(Fields[2]),std::stod(Fields[3]))); // C1
            if (count==3)   El.push_back(Vector3(std::stod(Fields[1]),std::stod(Fields[2]),std::stod(Fields[3]))); // C2
            if (count==4){
                            El.push_back(Vector3(std::stod(Fields[1]),std::stod(Fields[2]),std::stod(Fields[3]))); // C3
                            STL_Els_Raw.push_back(El);
                            El.clear();
            }
            count++;
            countglob++;
        }
        file.close();
        return true;
    }
    else{
        std::cout << FilePath << " . A suitable .stl file was not found at this position.\n";
        return false;
    }

}

//--- Geometry generation

void STL_Geometry::Generate_Nodes()
{
    // This temporarily stores the nodes of the STL import
    // Although almost all of these overlap, the STL lacks any connectivity information. So the all nodes must be created uniquely

    Inertial_CS = new CoordSys();

    Real Scale = 100.0;

    int IDCount = 0;
    for (int i=0; i<STL_Els_Raw.size(); i++){
        std::vector<SP_Node> STLNd;
        for (int j=0; j<STL_Els_Raw[i].size(); j++){
            Vector3 PScale = Scale*STL_Els_Raw[i][j];
            SP_Node N = std::make_shared<Node>(Inertial_CS,PScale);
            Nodes.push_back(N);
            STLNd.push_back(N);
            N->ID = IDCount;
            IDCount++;
        }
        STL_Els.push_back(STLNd);
    }

    STL_Els_Raw.clear();        // Clear list
}

void STL_Geometry::Generate_Elements()
{
    // This generates panel elements from the STL nodes

    for (int i=0; i<STL_Els.size(); i++)    Elements.push_back(std::make_shared<Tri_Element>(   STL_Els[i][2],
                                                                                                STL_Els[i][1],
                                                                                                STL_Els[i][0]));

    STL_Els.clear();            // Clear list to avoid problems with destructor
}

}
