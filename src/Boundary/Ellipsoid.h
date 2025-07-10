/****************************************************************************
    BEMUse Ellipsoid
    Copyright (C) 2022 Joseph Saverin j.saverin@tu-berlin.de

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Information on file:

    -> This generates an ellipsoid volume

*****************************************************************************/

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "Boundary_Base.h"
#include "Numerical_Wave_Tank.h"

namespace BEMUse
{

class Ellipsoid : public Boundary
{
protected:

    // Discretisation
    int NA, NZ;
    int NRFS=0, NAFS=0;     // Discretisation of free surface panels

    int NWall = 10;         // Discreitisation of wall panels

    // Geo vars
    Real a,b,c,Depth;
    Real RadFS = 0;         // Free-surface grid radius
    Real Floor = 0;         // Depth of sea bed / floor

public:

    //--- Constructor
    Ellipsoid()   {}

    //--- Geometry generation
    void Generate_Nodes() override;
    void Generate_Elements() override;
    void Generate_Aux_Nodes() override                               {}
    void Generate_Aux_Elements()    override                         {}

    //--- Geometry specification
    void Set_Parameters(std::vector<Parameter> &Params) override;

    //--- ID retrieval
    int Node_ID(int A, int Z);
};

class Semi_Ellipsoid : public Boundary      // Eg. Hemisphere
{

protected:

    // Discretisation
    int NA, NZ;
    int NRFS = 0;
    int NRES=0, NAES=0;
    int NWall = 10;
    // Geo vars
    Real a,b,c, RFS, Floor;

public:

    //--- Constructor
    Semi_Ellipsoid()   {}

    //--- Generate elements
    void Generate_Nodes();
    void Generate_Elements();
    void Generate_Aux_Nodes();
    void Generate_Aux_Elements();
    void Generate_FreeSurface_Nodes();
    void Generate_FreeSurface_Elements();

    //--- Geometry specification
    void Set_Parameters(std::vector<Parameter> &Params) override;

    //--- ID retrieval
    int Node_ID(int A, int Z);
    int Ext_Node_ID(int A, int Z);
    //    int Aux_Node_ID(int A, int Z)       {}
};

//--- Hybrid geometries: Ellipsoid submerged in a NWT

class Ellipsoid_NWT : public Boundary
{
protected:
    Boundary *Geo_Ellipsoid = nullptr;
    Boundary *Geo_Tank = nullptr;

public:
    //--- Constructor
    Ellipsoid_NWT()   {
        Geo_Ellipsoid = new Ellipsoid();
        Geo_Tank = new Numerical_Wave_Tank_Cylindrical();
    }

    // Geo setup
    void Set_Parameters(std::vector<Parameter> &Params) override {
        Geo_Ellipsoid->Set_Parameters(Params);
        Geo_Tank->Set_Parameters(Params);
    }

    void Setup() override {

        //--- Geometry is generated individually for both geos
        Geo_Ellipsoid->Generate_Nodes();
        Geo_Tank->Generate_Wall_Nodes();
        Geo_Tank->Generate_Floor_Nodes();
        Geo_Tank->Generate_FreeSurface_Nodes();

        Geo_Ellipsoid->Generate_Elements();
        Geo_Tank->Generate_Wall_Elements();
        Geo_Tank->Generate_Floor_Elements();
        Geo_Tank->Generate_FreeSurface_Elements();

        //--- Now absorb these into this object
        Geo_Ellipsoid->Get_Nodes(Nodes);
        Geo_Tank->Get_Wall_Nodes(Wall_Nodes);
        Geo_Tank->Get_Floor_Nodes(Floor_Nodes);
        Geo_Tank->Get_FreeSurface_Nodes(FreeSurface_Nodes);

        Geo_Ellipsoid->Get_Elements(Elements);
        Geo_Tank->Get_Wall_Elements(Wall_Elements);
        Geo_Tank->Get_Floor_Elements(Floor_Elements);
        Geo_Tank->Get_FreeSurface_Elements(FreeSurface_Elements);

        //--- Calculate geometric & kinematic quantities
        // for (SP_Geo G : Elements)               G->Set_Centroid();
        // for (SP_Geo G : Wall_Elements)          G->Set_Centroid();
        // for (SP_Geo G : Floor_Elements)         G->Set_Centroid();
        // for (SP_Geo G : FreeSurface_Elements)   G->Set_Centroid();

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

    // Destructor
    ~ Ellipsoid_NWT()
    {
        delete Geo_Ellipsoid;
        delete Geo_Tank;
    }
};

}

#endif // ELLIPSOID_H
