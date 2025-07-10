//---------------------------------------------------------------
//------------------ Ellipsoid Functions-------------------------
//---------------------------------------------------------------

#include "Ellipsoid.h"

namespace BEMUse
{

//--- Standard Ellipsoid

void Ellipsoid::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Semiaxis_a"))           a = P.Get_Param<Real>();
        if (P.myNameis("Semiaxis_b"))           b = P.Get_Param<Real>();
        if (P.myNameis("Semiaxis_c"))           c = P.Get_Param<Real>();
        if (P.myNameis("Ellipsoid_Depth"))      Depth = -P.Get_Param<Real>();
        if (P.myNameis("NPanels_Axial"))        NZ = P.Get_Param<int>();
        if (P.myNameis("NPanels_Azimuthal"))    NA = P.Get_Param<int>();
        if (P.myNameis("Cosine_Disc"))          Cosine = P.Get_Param<bool>();
        if (P.myNameis("Triangular_Panels"))    TriPanels = P.Get_Param<bool>();
    }

    // if (TriPanels) Panels_Outward = true;       // Hack for Rankine source
}

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
            if (TriPanels) {
                if (i%2==0) u += 0.5*TwoPI/NAR;
            }
            //------------------------------

            Vector3 BP  ( a*cos(u)*sinv, b*sin(u)*sinv, c*cosv);
            Vector3 dPdu(-a*sin(u)*sinv, b*cos(u)*sinv, 0);         dPdu.normalize();
            // Vector3 dPdv( a*cos(u)*cosv, b*sin(u)*cosv,-c*sinv);    dPdv.normalize();    // Normals pointing outward
            Vector3 dPdv(-a*cos(u)*cosv,-b*sin(u)*cosv, c*sinv);    dPdv.normalize();       // Normals pointing inward
            // Vector3 dPdv(-a*cos(u)*cosv,-b*sin(u)*cosv,-c*sinv);    dPdv.normalize();       // Normals pointing inward
            if (i==0){
                dPdu = Vector3(0,1,0);
                // dPdv = Vector3(-1,0,0);      // Normals pointing outward
                dPdv = Vector3(1,0,0);          // Normals pointing inward
            }
            if (i==NZ){
                dPdu = Vector3(0,1,0);
                // dPdv = Vector3(1,0,0);       // Normals pointing outward
                dPdv = Vector3(-1,0,0);         // Normals pointing inward
            }
            Vector3 dPdn = dPdu.cross(dPdv);                        dPdn.normalize();
            Quat tO = Quat::FromTwoVectors(UnitZ,dPdn);

            Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,BP));
        }
    }

    // Set node ID
    for (int i=0; i<Nodes.size(); i++)  Nodes[i]->ID = i;

    std::cout << "NNodes = " << size(Nodes) << std::endl;
}

void Ellipsoid::Generate_Elements()
{
    // This is the function to generate the surface elements over the geometry
    // These are created before creatin panels based on them.

    for (int i=0; i<NA; i++){
        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                         Nodes[Node_ID(i+1,1)],
                                                         Nodes[Node_ID(i,1)]));
        if (Panels_Outward) Elements.back()->Reorder_Nodes();
    }

    if (TriPanels) {
        // Generate central elements with triangular panels
        for (int z=1; z<NZ-1; z++){
            for (int i=0; i<NA; i++){
                if (z%2==0){
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i,z)],
                                                                     Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i+1,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i+2,z+1)],
                                                                     Nodes[Node_ID(i+1,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                }
                else
                {
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i,z)],
                                                                     Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i+1,z+1)],
                                                                     Nodes[Node_ID(i,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                }
            }
        }
    }
    else
    {
        // Generate central elements with quadratic panels
        for (int z=1; z<NZ-1; z++){
            // for (int z=1; z<NZ/2; z++){
            for (int i=0; i<NA; i++){
                Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,z)],
                                                                  Nodes[Node_ID(i+1,z)],
                                                                  Nodes[Node_ID(i+1,z+1)],
                                                                  Nodes[Node_ID(i,z+1)]));
                if (Panels_Outward) Elements.back()->Reorder_Nodes();
            }
        }
    }

    // Nose (triangular elements)
    for (int i=0; i<NA; i++){
        Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,NZ-1)],
                                                         Nodes[Node_ID(i,NZ)],
                                                         Nodes[Node_ID(i,NZ-1)]));
        if (Panels_Outward) Elements.back()->Reorder_Nodes();
    }

    for (int i=0; i<Elements.size(); i++)  Elements[i]->Set_Centroid();
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

void Semi_Ellipsoid::Set_Parameters(std::vector<Parameter> &Params)
{
    // This sets the parameters for the geometry and simulation
    StdAppend(Parameters, Params);
    for (Parameter P : Parameters)
    {
        if (P.myNameis("Semiaxis_a"))                       a = P.Get_Param<Real>();
        if (P.myNameis("Semiaxis_b"))                       b = P.Get_Param<Real>();
        if (P.myNameis("Semiaxis_c"))                       c = P.Get_Param<Real>();
        if (P.myNameis("FreeSurface_Radius"))               RFS = P.Get_Param<Real>();
        if (P.myNameis("NPanels_Axial"))                    NZ = P.Get_Param<int>();
        if (P.myNameis("NPanels_Azimuthal"))                NA = P.Get_Param<int>();
        if (P.myNameis("NPanels_FreeSurface_Radial"))       NRES = P.Get_Param<int>();
        if (P.myNameis("NPanels_FreeSurface_Int_Radial"))   NRFS = P.Get_Param<int>();
        if (P.myNameis("Cosine_Disc"))                      Cosine = P.Get_Param<bool>();
        if (P.myNameis("Triangular_Panels"))                TriPanels = P.Get_Param<bool>();
    }

    // if (TriPanels) Panels_Outward = true;       // Hack for Rankine source

    // Some extra parameters for the Rankine source test case
    if (RFS<=1.0)   RFS = 2.;       // Ensure FS radius is sufficiently large
    NAES = NA;                      // Ensure the azimuthal discretisation of the FS is the same as the ellipsoid
}

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

    if (TriPanels) NA = NAES;

    //---- Generate coordinates

    for (int i=0; i<NZ+1; i++)          // Zenith angle v must sweep from PI-PI/2
    {
        Real v;
        if (Cosine)    v = PI*(1.0-0.5*PiCosFac(i,NZ+1));
        else           v = PI - 0.5*PI*i/NZ;
        Real sinv = sin(v), cosv = cos(v);
        if (i==0)       { sinv = 0.0;   cosv = -1.0;}
        else if (i==NZ) { sinv = 1.0;   cosv = 0.0;}

        int NAR = NA;
        if (i==0)               NAR = 1;
        for (int j = 0; j<NAR; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NAR;

            //------Mod for triangular panels---------
            if (TriPanels) {
                if (i%2==0) u += 0.5*TwoPI/NAR;
            }
            //----------------------------------------
            Vector3 BP  ( a*cos(u)*sinv, b*sin(u)*sinv, c*cosv);
            Vector3 dPdu(-a*sin(u)*sinv, b*cos(u)*sinv, 0);         dPdu.normalize();
            // Vector3 dPdv( a*cos(u)*cosv, b*sin(u)*cosv,-c*sinv);    dPdv.normalize();    // Normals pointing outward
            Vector3 dPdv(-a*cos(u)*cosv,-b*sin(u)*cosv, c*sinv);    dPdv.normalize();       // Normals pointing inward

            if (i==0){
                dPdu = Vector3(0,1,0);
                // dPdv = Vector3(-1,0,0);
                dPdv = Vector3(1,0,0);          // Normals pointing inward
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
    for (int i=0; i<NA; i++){
        Elements.push_back(std::make_shared<Tri_Element>(   Nodes[Node_ID(0,0)],
                                                         Nodes[Node_ID(i+1,1)],
                                                         Nodes[Node_ID(i,1)]));
        if (Panels_Outward) Elements.back()->Reorder_Nodes();
    }

    if (TriPanels) {
        // Generate central elements with triangular panels
        for (int z=1; z<NZ; z++){
            for (int i=0; i<NA; i++){
                if (z%2==0){
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i,z)],
                                                                     Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i+1,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i+2,z+1)],
                                                                     Nodes[Node_ID(i+1,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                }
                else
                {
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i,z)],
                                                                     Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                    Elements.push_back(std::make_shared<Tri_Element>(Nodes[Node_ID(i+1,z)],
                                                                     Nodes[Node_ID(i+1,z+1)],
                                                                     Nodes[Node_ID(i,z+1)]));
                    if (Panels_Outward) Elements.back()->Reorder_Nodes();
                }
            }
        }
    }
    else
    {
        // Generate quadratic elements
        for (int z=1; z<NZ; z++){
            for (int i=0; i<NA; i++){
                Elements.push_back(std::make_shared<Quad_Element>(  Nodes[Node_ID(i,z)],
                                                                  Nodes[Node_ID(i+1,z)],
                                                                  Nodes[Node_ID(i+1,z+1)],
                                                                  Nodes[Node_ID(i,z+1)]));
                if (Panels_Outward) Elements.back()->Reorder_Nodes();
            }
        }
    }

    for (int i=0; i<Elements.size(); i++)  Elements[i]->Set_Centroid();
}

void Semi_Ellipsoid::Generate_Aux_Nodes()
{
    // This creates the nodes at the free surface

    if (TriPanels) return;      // Jump out: Avoid for Rankine Source test case.

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

    if (TriPanels) return;      // Jump out: Avoid for Rankine Source test case.

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

void Semi_Ellipsoid::Generate_FreeSurface_Nodes()
{
    // A volume of revolution has as it's surface a circle. We simply generate an additional array of elements

    // Factor...
    // The factor RFS scales the fre surfae based on the dimensions of the ellipsoid
    // The scaling factor goes from 1 - RFS. and varies either linearly or cosine to increase resolution near the hemisphere.
    Real FShift = 0.;
    if (!TriPanels) FShift = 0.001; // This is added in to avoid issues with FS elevation

    for (int i=0; i<NRES+1; i++)
    {
        Real f;
        if (Cosine) f = 1.0 + FShift + (RFS-1.0)*PiCosFac(i,NRES+1);     // Cosine
        else        f = 1.0 + FShift + (RFS-1.0)*i*1.0/NRES;             // Linear

        for (int j = 0; j<NA; j++)     // Azimuth angle u must sweep from 0-2PI
        {
            Real u = j*TwoPI/NA;

            //--- Mod for triangular panels
            if (TriPanels) {if (i%2==0) u += 0.5*TwoPI/NA;}
            Vector3 BP  ( a*f*cos(u), b*f*sin(u), 0);
            Quat tO = Quat::FromTwoVectors(UnitZ,Vector3(0.,0.,-1.));
            FreeSurface_Nodes.push_back(std::make_shared<Node>(Inertial_CS,tO,BP));
        }
    }

    for (int i=0; i<FreeSurface_Nodes.size(); i++) FreeSurface_Nodes[i]->ID = i;
}

void Semi_Ellipsoid::Generate_FreeSurface_Elements()
{
    // Generate the panels on the free surface exterior to the body

    if (TriPanels){

        for (int z=0; z<NRES; z++){
            for (int i=0; i<NA; i++){
                if (z%2==0){
                    FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i,z)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i+1,z)]));
                    FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i+1,z)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i+2,z+1)]));
                }
                else
                {
                    FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i,z)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i,z+1)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i+1,z)]));
                    FreeSurface_Elements.push_back(std::make_shared<Tri_Element>(FreeSurface_Nodes[Ext_Node_ID(i+1,z)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i,z+1)],
                                                                                 FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)]));
                }
            }
        }
    }
    else {

        for (int z=0; z<NRES; z++){
            for (int i=0; i<NA; i++)  FreeSurface_Elements.push_back(std::make_shared<Quad_Element>(FreeSurface_Nodes[Ext_Node_ID(i,z+1)],
                                                                              FreeSurface_Nodes[Ext_Node_ID(i+1,z+1)],
                                                                              FreeSurface_Nodes[Ext_Node_ID(i+1,z)],
                                                                              FreeSurface_Nodes[Ext_Node_ID(i,z)]));
        }

    }

    for (int i=0; i<FreeSurface_Elements.size(); i++)  FreeSurface_Elements[i]->Set_Centroid();
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
