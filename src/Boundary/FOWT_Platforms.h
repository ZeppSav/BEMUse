/****************************************************************************
    BEMUse FOWT Platforms
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

    -> This contains the classes which specifically define offshore wind turbine platforms

*****************************************************************************/

#ifndef FOWT_PLATFORMS_H
#define FOWT_PLATFORMS_H

#include "Volume_of_Revolution.h"

namespace BEMUse
{

class OC3_SparBuoy : public Tapered_SparBuoy
{
protected:

    //--- Geo parameters
    Real ScaleFac = 1.0;

public:

    //--- Constructor
    OC3_SparBuoy()  {}

    //--- Geometry specification
    void Set_Parameters(std::vector<Parameter> &Params) override
    {
        StdAppend(Parameters, Params);
        for (Parameter P : Parameters)
        {
            if (P.myNameis("Scaling_Factor"))                   ScaleFac = P.Get_Param<Real>();
            if (P.myNameis("FreeSurface_Radius"))               RFS = P.Get_Param<Real>();
            if (P.myNameis("NPanels_Radial"))                   NR = P.Get_Param<int>();
            if (P.myNameis("NPanels_Vertical_Section1"))        NV1 = P.Get_Param<int>();
            if (P.myNameis("NPanels_Vertical_Section2"))        NV2 = P.Get_Param<int>();
            if (P.myNameis("NPanels_Vertical_Section3"))        NV3 = P.Get_Param<int>();
            if (P.myNameis("NPanels_Azimuthal"))                NA = P.Get_Param<int>();
            if (P.myNameis("NPanels_FreeSurface_Radial"))       NRES = P.Get_Param<int>();
            if (P.myNameis("NPanels_FreeSurface_Int_Radial"))   NRFS = P.Get_Param<int>();
            if (P.myNameis("Cosine_Disc"))                      Cosine = P.Get_Param<bool>();
            if (P.myNameis("Triangular_Panels"))                TriPanels = P.Get_Param<bool>();
        }

        // Set floater dimensions
        RB = 4.7*ScaleFac;
        RT = 3.25*ScaleFac;
        Z = -120*ScaleFac;
        H1 = 108*ScaleFac;
        H2 = 8*ScaleFac;
        R = RT;
    }

    //--- Geometry specification
    void Set_Dimensions(std::vector<Real> &D)
    {
        ScaleFac = D[0];
        RB = 4.7*ScaleFac;
        RT = 3.25*ScaleFac;
        Z = -120*ScaleFac;
        H1 = 108*ScaleFac;
        H2 = 8*ScaleFac;
    }

    void Set_External_Dimensions(std::vector<Real> &D)      {RFS = D[0]; R = RT;}

    //--- Geometry calculations (Hard coded)
    void Set_Kin_Params()
    {
        // Hard-code here the kinematic parameters
        Centre_Gravity = Vector3(0,0,-89.9155);
        Mass = 7.46633e+06;

        Mass_Mat <<         7.46633e+06,  	0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   	7.46633e+06,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   	0.00000e+00,   7.46633e+06,   0.00000e+00,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   	0.00000e+00,   0.00000e+00,   4.22923e+09,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   	0.00000e+00,   0.00000e+00,   0.00000e+00,   4.22923e+09,   0.00000e+00,
                            0.00000e+00,   	0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   1.226e+10;

        Stiffness_Mat <<    0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   0.00000e+00,   3.32941e+05,    0.00000e+00,  0.00000e+00,   0.00000e+00,
                            0.00000e+00,   0.00000e+00,   0.00000e+00,  -4.99918e+09,   0.00000e+00,   0.00000e+00,
                            0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,  -4.99918e+09,   0.00000e+00,
                            0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00,   0.00000e+00;

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

};

class Softwind_Upscaled : public Tapered_SparBuoy
{
protected:

public:
    //--- Constructor
    Softwind_Upscaled() {}

    //--- Geometry specification
    void Set_Dimensions(std::vector<Real> &D)
    {
        RB = 9.0;
        RT = 5.6;
        Z = -91.4;
        H1 = 78;
        H2 = 8;
        RFS = RT;
    }

    //--- Geometry calculations (Hard coded)
    void Set_Kin_Params()
    {
        // Hard-code here the kinematic parameters
        // Scale softwind lab to experiment:    Scale 1:40
        // Keel Depth: 2.285   (exp)   91.4 (full scale)
        //
        // Centre of gravity (taken from keel): 0.496 -> Full scale = -91.4 + 40*0.496 = -71.56
        Mass = 1.94E+07;
        Centre_Gravity = Vector3(0, 0, -71.56);
        // Need to specify radii of gyration
    }
};

}

#endif // FOWT_PLATFORMS_H
