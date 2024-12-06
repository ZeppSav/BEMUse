//---------------------------------------------------------------
//------------------ Ellipsoid Functions-------------------------
//---------------------------------------------------------------

#include "Ellipsoid.h"

namespace BEMUse
{

//--- Standard Ellipsoid

void Ellipsoid::Generate_Nodes()
{

    // Create coordinates of an ellipsoid
    // The surface of an ellipsoid of side a,b.
    // r =      [ a cos(u)sin(v), b sin(u)sin(v), c cos v];       0 < u < 2*PI, PI/2 < z < PI
    // dr/du =  [-a sin(u)sin(v), b cos(u)sin(v), c cos v];       0 < u < 2*PI, PI/2 < z < PI
    // dr/dv =  [ a cos(u)cos(v), b sin(u)cos(v),-c sin v];       0 < u < 2*PI, PI/2 < z < PI
    // NA, NZ represents the discretisation of the surface
    // The origin of the coordinate system is at [X,Y,Z] = [0,0,0]

    if (!Global_CS)    Global_CS = new CoordSys();
    //    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);

    Vector3 P(0,0,Depth);
    Quat O(Quat::Identity());
    //    Quat O(Eigen::AngleAxisf(0.5*PI, UnitY));// = Quat(Eigen::AngleAxisf(0.5*PI, UnitY));
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS,O,P); // Rotated

    //---- Generate coordinates

    for (int i=0; i<NZ+1; i++)          // Zenith angle v must sweep from PI-PI/2
    {
        Real v;
        if (Cosine) v = PI*(1.0-1.0*PiCosFac(i,NZ+1));   // Half cosine dist
        else        v = PI - PI*i/(NZ);                    // Linear cosine dist (sphere)
        Real sinv = sin(v), cosv = cos(v);
        if (i==0)       { sinv = 0.0;   cosv = -1.0;}
        else if (i==NZ) { sinv = 0.0;   cosv = 1.0;}

        int NAR = NA;
        if (i==0 || i==NZ)      NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;

            //--- Mod for triangular panels
            if (i%2==0) u += 0.5*TwoPI/NAR;
            //------------------------------

            Vector3 BP  ( a*cos(u)*sinv, b*sin(u)*sinv, c*cosv);
            Vector3 dPdu(-a*sin(u)*sinv, b*cos(u)*sinv, 0);         dPdu.normalize();
            Vector3 dPdv( a*cos(u)*cosv, b*sin(u)*cosv,-c*sinv);    dPdv.normalize();
            if (i==0){
                dPdu = Vector3(0,1,0);
                dPdv = Vector3(-1,0,0);
            }
            if (i==NZ){
                dPdu = Vector3(0,1,0);
                dPdv = Vector3(1,0,0);
            }
            Vector3 dPdn = dPdu.cross(dPdv);                        dPdn.normalize();
            Quat tO = Quat::FromTwoVectors(UnitZ,dPdn);

            Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,BP));
        }
    }

    // Set node ID
    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

    //    // Calculate analytical solution parameters
    //    int Xmax = 100;
    //    Real dX = 0.00001;
    //    alpha0 = 0.0;
    //    beta0 = 0.0;
    //    for (Real x=0.5*dX; x<Xmax; x+=dX)
    //    {
    //        Real delconst = a*b*c*dX/sqrt((a*a+x)*(b*b+x)*(c*c+x));
    //        alpha0 += delconst/(a*a+x);
    //        beta0 += delconst/(b*b+x);
    //    }
    //    gamma0 = 2.0-alpha0-beta0;
    //    std::cout << "Ellipsoid constants given by: alpha0 = " << alpha0 << ", beta0 = " << beta0 << ",  gamma0 = " << gamma0 << std::endl;
}

void Ellipsoid::Generate_Elements()
{
    // This is the function to generate the surface elements over the geometry
    // These are created before creatin panels based on them.

    std::cout << "Generating Surface Elements \n";

    // Base (triangular elements)
    for (int i=0; i<NA; i++)        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                                                        Nodes[Node_ID(i+1,1)],
                                                                                        Nodes[Node_ID(i,1)]));

    // Central quadratic elements
    // for (int z=1; z<NZ-1; z++){
    //     for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,z)],
    //                                                                                     Nodes[Node_ID(i+1,z)],
    //                                                                                     Nodes[Node_ID(i+1,z+1)],
    //                                                                                     Nodes[Node_ID(i,z+1)]));
    // }

    //--- Mod for triangular panels
    for (int z=1; z<NZ-1; z++){
        for (int i=0; i<NA; i++){
             if (z%2==0){
                Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i,z)],
                                                                 Nodes[Node_ID(i+1,z+1)],
                                                                 Nodes[Node_ID(i+1,z)]));
                Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,z)],
                                                                 Nodes[Node_ID(i+1,z+1)],
                                                                 Nodes[Node_ID(i+2,z+1)]));
             }
             else
             {
                 Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i,z)],
                                                                  Nodes[Node_ID(i,z+1)],
                                                                  Nodes[Node_ID(i+1,z)]));
                 Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,z)],
                                                                  Nodes[Node_ID(i,z+1)],
                                                                  Nodes[Node_ID(i+1,z+1)]));
             }
        }
    }
    //------------------------------

    // Nose (triangular elements)
    for (int i=0; i<NA; i++)        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(i+1,NZ-1)],
                                                                                        Nodes[Node_ID(i,NZ)],
                                                                                        Nodes[Node_ID(i,NZ-1)]));

//    for (int i=0; i<NA; i++)        std::cout <<Node_ID(i+1,NZ-1) << " " << Node_ID(i,NZ) << " "<< Node_ID(i,NZ-1) << std::endl;

//        std::cout << "Size: " <<  Nodes.size() <<" "<< Elements.size() << std::endl;
}

int Ellipsoid::Node_ID(int A, int Z)
{
    // This returns the node ID number for the given azimuthal and axial position

    if (Z==NZ)          return 1+NA*(NZ-1);
    if (Z==0)           return 0;
    else {
        if (A>=NA)      return 1+(Z-1)*NA+(A-NA);
//        else if (A<0)   return 1+(Z-1)*NA+(A+NA);
        else            return 1+(Z-1)*NA+A;
    }
}

//--- Semi Ellipsoid

void Semi_Ellipsoid::Generate_Nodes()
{

    // Create coordinates of an ellipsoid
    // The surface of an ellipsoid of side a,b.
    // r =      [ a cos(u)sin(v), b sin(u)sin(v), c cos v];       0 < u < 2*PI, PI/2 < z < PI
    // dr/du =  [-a sin(u)sin(v), b cos(u)sin(v), c cos v];       0 < u < 2*PI, PI/2 < z < PI
    // dr/dv =  [ a cos(u)cos(v), b sin(u)cos(v),-c sin v];       0 < u < 2*PI, PI/2 < z < PI
    // NA, NZ represents the discretisation of the surface
    // The origin of the coordinate system is at [X,Y,Z] = [0,0,0]

    if (!Global_CS)    Global_CS = new CoordSys();
    if (!Inertial_CS)  Inertial_CS = new CoordSys(Global_CS);
    //---- Generate coordinates

    for (int i=0; i<NZ+1; i++)          // Zenith angle v must sweep from PI-PI/2
    {
        Real v;
        if (Cosine)    v = PI*(1.0-0.5*PiCosFac(i,NZ+1));
        else           v = PI - 0.5*PI*i/NZ;
//        Real v = PI*(1.0-0.5*PiCosFac(i,NZ+1));
//            Real v = PI - 0.5*PI*i/(NZ);                    // Linear cosine dist
        Real sinv = sin(v), cosv = cos(v);
        if (i==0)       { sinv = 0.0;   cosv = -1.0;}
        else if (i==NZ) { sinv = 1.0;   cosv = 0.0;}

        int NAR = NA;
        if (i==0)               NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;
            Vector3 BP  ( a*cos(u)*sinv, b*sin(u)*sinv, c*cosv);
            Vector3 dPdu(-a*sin(u)*sinv, b*cos(u)*sinv, 0);         dPdu.normalize();
            Vector3 dPdv( a*cos(u)*cosv, b*sin(u)*cosv,-c*sinv);    dPdv.normalize();
            if (i==0){
                dPdu = Vector3(0,1,0);
                dPdv = Vector3(-1,0,0);
            }
            Vector3 dPdn = dPdu.cross(dPdv);                        dPdn.normalize();
            Quat tO = Quat::FromTwoVectors(UnitZ,dPdn);
            Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,BP));
        }
    }

    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;
}

void Semi_Ellipsoid::Generate_Elements()
{
    // This is the function to generate the surface elements over the geometry
    // These are created before creatin panels based on them.

    std::cout << "Generating Surface Elements \n";

    // Base (triangular elements)
    for (int i=0; i<NA; i++)        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                                                        Nodes[Node_ID(i+1,1)],
                                                                                        Nodes[Node_ID(i,1)]));

    // Central quadratic elements
    for (int z=1; z<NZ; z++){
        for (int i=0; i<NA; i++)    Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,z)],
                                                                                        Nodes[Node_ID(i+1,z)],
                                                                                        Nodes[Node_ID(i+1,z+1)],
                                                                                        Nodes[Node_ID(i,z+1)]));
    }
}

void Semi_Ellipsoid::Generate_Aux_Nodes()
{
    // This creates the nodes at the free surface
    for (int i=0; i<NRFS+1; i++)          // Zenith angle v must sweep from PI-PI/2
    {
        Real f;
        if (Cosine)    f = PiCosFac(i,NRFS+1);   // Cosine
        else           f = 1.0*i/NRFS;

        int NAR = NA;
        if (i==0)      NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;
            Vector3 BP  ( a*f*cos(u), b*f*sin(u), 0);
            Aux_Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    for (int i=0; i<Aux_Nodes.size(); i++)  Aux_Nodes[i]->ID = i;
}

void Semi_Ellipsoid::Generate_Aux_Elements()
{
    // This is the function to generate the surface elements over the geometry
    // These are created before creating panels based on them.

    // Base (triangular elements)
    for (int i=0; i<NA; i++)        Aux_Elements.push_back(std::make_shared<Tri_Element>(   Aux_Nodes[Node_ID(0,0)],
                                                                                            Aux_Nodes[Node_ID(i,1)],
                                                                                            Aux_Nodes[Node_ID(i+1,1)]));

    // Central quadratic elements
    for (int z=1; z<NRFS; z++){
        for (int i=0; i<NA; i++)    Aux_Elements.push_back(std::make_shared<Quad_Element>(  Aux_Nodes[Node_ID(i,z+1)],
                                                                                            Aux_Nodes[Node_ID(i+1,z+1)],
                                                                                            Aux_Nodes[Node_ID(i+1,z)],
                                                                                            Aux_Nodes[Node_ID(i,z)]));
    }

//        for (SP_Geo G : Aux_Elements)   G->Set_Centroid();
}

void Semi_Ellipsoid::Generate_Ext_Nodes()
{
    // A volume of revolution has as it's surface a circle. We simply generate an additional array of elements

//    NRES = 2*NRFS;  // Hack!
//    NAES = 2*NA;

    for (int i=0; i<NRES+1; i++)
    {
        Real f = (1.001+RFS*i/NRES);           // Linear
//        Real r = RFS*PiCosFac(i,NRES+1);   // Cosine
        for (int j = 0; j<NAES; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAES;
//            Vector3 BP  ( a*f*cos(u) + Origin(0), b*f*sin(u) + Origin(1), Origin(2));
            Vector3 BP  ( a*f*cos(u), b*f*sin(u), 0);
            Ext_Nodes.push_back(std::make_shared<Node>(Inertial_CS,BP));
        }
    }

    for (int i=0; i<Ext_Nodes.size(); i++) Ext_Nodes[i]->ID = i;
}

void Semi_Ellipsoid::Generate_Ext_Elements()
{
    // Central quadratic elements
    for (int z=0; z<NRES; z++){
        for (int i=0; i<NAES; i++)    Ext_Elements.push_back(std::make_shared<Quad_Element>(    Ext_Nodes[Ext_Node_ID(i,z+1)],
                                                                                                Ext_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                                                Ext_Nodes[Ext_Node_ID(i+1,z)],
                                                                                                Ext_Nodes[Ext_Node_ID(i,z)]));
    }

    for (int i=0; i<Ext_Elements.size(); i++)  Ext_Elements[i]->Set_Centroid();
}

int Semi_Ellipsoid::Node_ID(int A, int Z)
{
    // This returns the node ID number for the given azimuthal and axial position

    if (Z==0) return 0;
    else {
        if (A>=NA)      return 1+(Z-1)*NA+(A-NA);
        else if (A<0)   return 1+(Z-1)*NA+(A+NA);
        else            return 1+(Z-1)*NA+A;
    }
}

int Semi_Ellipsoid::Ext_Node_ID(int A, int Z)
{
    // Returns the node ID based on the azimuthal and axial position
    if (A>=NAES)    return Z*NAES + (A-NAES);
    else            return Z*NAES + A;
}

}
