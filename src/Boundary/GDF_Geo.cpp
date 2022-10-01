//-----------------------------------------------------------------------------
//-------------------------GDF Geo routines------------------------------------
//-----------------------------------------------------------------------------

#include "GDF_Geo.h"

namespace BEMUse
{

//--- Import/Export

bool GDF_Geometry::Read_Input_File(std::string &FilePath)
{
    // This function opens the File with filename in the compile directory and reads
    // through the file storing node positions and panel connectivity

    GDF_Els_Raw.clear();

    std::ifstream file(FilePath);
    std::string line;
    if (file.is_open())
    {
        int count=0;
        std::vector<Vector3> El;

        while ( std::getline (file,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields.size()==3)
            {
                // Read in elements
                if (count==0)   El.push_back(Vector3(std::stod(Fields[0]),std::stod(Fields[1]),std::stod(Fields[2])));
                if (count==1)   El.push_back(Vector3(std::stod(Fields[0]),std::stod(Fields[1]),std::stod(Fields[2])));
                if (count==2)   El.push_back(Vector3(std::stod(Fields[0]),std::stod(Fields[1]),std::stod(Fields[2])));
                if (count==3)
                {
                                Vector3 P4(std::stod(Fields[0]),std::stod(Fields[1]),std::stod(Fields[2]));
                                if (P4 != El[2])    El.push_back(P4);      // If quadrilateral, add. Otherwise ignore.
                                GDF_Els_Raw.push_back(El);
                                El.clear();
                }

                // Increment count
                if (count==3)   count=0;
                else            count++;
             }

        }
        file.close();
        return true;
    }
    else{
        std::cout << FilePath << " . A suitable .gdf file was not found at this position.\n";
        return false;
    }

}

//--- Geometry generation

void GDF_Geometry::Generate_Nodes()
{
    // This temporarily stores the nodes of the GDF import
    // Although almost all of these overlap, the GDF lacks any connectivity information. So the all nodes must be created uniquely

    Inertial_CS = new CoordSys();

    Real Scale = 1.0;

    int IDCount = 0;
    for (int i=0; i<GDF_Els_Raw.size(); i++){
        std::vector<SP_Node> STLNd;
        for (int j=0; j<GDF_Els_Raw[i].size(); j++){
            Vector3 PScale = Scale*GDF_Els_Raw[i][j];
            SP_Node N = std::make_shared<Node>(Inertial_CS,PScale);
            Nodes.push_back(N);
            STLNd.push_back(N);
            N->ID = IDCount;
            IDCount++;
        }
        GDF_Els.push_back(STLNd);
    }

    GDF_Els_Raw.clear();        // Clear list
}

void GDF_Geometry::Generate_Elements()
{
    // This generates panel elements from the GDF nodes

    for (int i=0; i<GDF_Els.size(); i++){

        if (GDF_Els[i].size()==3)   Elements.push_back(std::make_shared<Tri_Element>(   GDF_Els[i][0],
                                                                                        GDF_Els[i][1],
                                                                                        GDF_Els[i][2]));
        if (GDF_Els[i].size()==4)   Elements.push_back(std::make_shared<Quad_Element>(  GDF_Els[i][0],
                                                                                        GDF_Els[i][1],
                                                                                        GDF_Els[i][2],
                                                                                        GDF_Els[i][3]));
    }

    GDF_Els.clear();            // Clear list to avoid problems with destructor
}

}
