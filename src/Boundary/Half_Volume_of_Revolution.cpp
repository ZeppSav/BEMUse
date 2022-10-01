//---------------------------------------------------------------
//------------------ HVOR Functions------------------------------
//---------------------------------------------------------------

#include "Half_Volume_of_Revolution.h"

namespace BEMUse
{

//--- Geometry functions

void Half_Volume_of_Revolution::Generate_Nodes()
{
    // This generates the ndoes for the semi cylinder

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    // Create perimeter
    Set_Perimeter();

    // End face 1
    for (int i=0; i<NR+1; i++){

        Real RF;
        if (Cosine) RF = PiCosFac(i,NR+1);   // Cosine;
        else        RF = 1.0*i/(NR+1);

        if (i==0)    Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(-0.5*L,0,0)));
        else{
            for (int j=0; j<NA+1; j++){
                Vector3 BP(-0.5*L,RF*Perimeter[j](1),RF*Perimeter[j](2));
                Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
            }
        }
    }

    // Create cylindrical section
    for (int i=0; i<NL+1; i++){

        Real tL;
        if (Cosine) tL = -0.5*L*(1.0-2.0*PiCosFac(i,NL+1));   // Cosine;
        else        tL = -0.5*L*(1.0-2.0*i/NL);

        for (int j=0; j<NA+1; j++){
            Vector3 BP(tL,Perimeter[j](1),Perimeter[j](2));
            Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    // End face 2
    for (int i=1; i<NR+1; i++){

        Real RF;
        if (Cosine) RF = R*(1.0-PiCosFac(i,NR+1));   // Cosine;
        else        RF = R*(1.0-1.0*i/(NR+1));

        if (i==NR)    Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(0.5*L,0,0)));
        else{
            for (int j=0; j<NA+1; j++){
                Vector3 BP(0.5*L,RF*Perimeter[j](1),RF*Perimeter[j](2));
                Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
            }
        }
    }

    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

//    for (int i=0; i<Nodes.size(); i++) {
//        Vector3 P = Nodes[i]->Position_Global();
//        std::cout << P(0) << " " << P(1) << " " << P(2) << std::endl;
//    }
 }

void Half_Volume_of_Revolution::Generate_Aux_Nodes()
{
    // Interior free surface points
    for (int i=0; i<NXAS+1; i++){
        for (int j=0; j<NYAS+1; j++){
            Real x,y;
            if (Cosine){
                x = -0.5*L + L*PiCosFac(i,NXAS+1);
                y = -R + 2.0*R*PiCosFac(j,NYAS+1);
            }
            else {
                x = -0.5*L + L*i/NXAS;
                y = -R + 2.0*R*j/NYAS;
            }
            Aux_Nodes.push_back(std::make_shared<Node>(Inertial_CS,Vector3(x,y,0.0)));
        }
    }

    for (int i=0; i<Aux_Nodes.size(); i++)  Aux_Nodes[i]->ID = i;
}

void Half_Volume_of_Revolution::Generate_Elements()
{
    // Creates panels

    // Face 1 panels
    for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                                                    Nodes[Node_ID(i,1)],
                                                                                    Nodes[Node_ID(i+1,1)]));

    // First the end panels
    for (int j=1; j<NR; j++){
        for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,j)],
                                                                                        Nodes[Node_ID(i,j+1)],
                                                                                        Nodes[Node_ID(i+1,j+1)],
                                                                                        Nodes[Node_ID(i+1,j)] ));
    }

    for (int j=NR; j<NR+NL; j++){
        for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,j)],
                                                                                        Nodes[Node_ID(i,j+1)],
                                                                                        Nodes[Node_ID(i+1,j+1)],
                                                                                        Nodes[Node_ID(i+1,j)] ));
    }

    for (int j=NR+NL; j<2*NR+NL-1; j++){
        for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,j)],
                                                                                        Nodes[Node_ID(i,j+1)],
                                                                                        Nodes[Node_ID(i+1,j+1)],
                                                                                        Nodes[Node_ID(i+1,j)] ));
    }

    // Face 2 panels
    for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,2*NR+NL)],
                                                                                    Nodes[Node_ID(i+1,2*NR+NL-1)],
                                                                                    Nodes[Node_ID(i,2*NR+NL-1)]));
}

void Half_Volume_of_Revolution::Generate_Aux_Elements()
{
    // First the end panels
    for (int i=0; i<NXAS; i++){
        for (int j=0; j<NYAS; j++)   Aux_Elements.push_back(std::make_shared<Quad_Element>( Aux_Nodes[(i+1)*(NYAS+1) + j],
                                                                                            Aux_Nodes[(i+1)*(NYAS+1) + j+1],
                                                                                            Aux_Nodes[i*(NYAS+1) + j+1],
                                                                                            Aux_Nodes[i*(NYAS+1) + j] ));
    }
}

void Half_Volume_of_Revolution::Set_Perimeter()
{
    // The default perimeter is a half cylinder

    for (int i=0; i<NA+1; i++){
        Real TH;
        if (Cosine) TH = PI*PiCosFac(i,NA+1);   // Cosine;
        else        TH = PI*i/(NA);
        Perimeter.push_back(Vector3(0.0,R*cos(TH),-R*sin(TH)));
    }
}

int Half_Volume_of_Revolution::Node_ID(int A, int Z)
{
    if (Z==0) return 0;
    else if (Z==(2*NR+NL)) return Nodes.size()-1;
    else return 1+Z*(NA+1)+A;
}

}
