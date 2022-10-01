/****************************************************************************
    BEMUse Library
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

    -> A set of function which provide the user with information about BEMUse

*****************************************************************************/

#include "BEMUse_IO.h"

#ifndef BEMUSE_INQUIRY_H
#define BEMUSE_INQUIRY_H

bool BEMUse_Enquiry(int argc, char *argv[])
{
    // This function essentially passes through all possible options of the console inputs
    // and responds appropriately.

    if (argc==1){
        std::cout << "\n";
        std::cout << "Good news! BEMUse is in your directory.\n";
        std::cout << "Bad news... BEMUse needs more inputs to be helpful.\n";
        return true;
    }

    // Information has been requested in some form: Convert these to strings
    std::string S1, S2, S3;

    if (argc>1) S1 = std::string(argv[1]);
    if (argc>2) S2 = std::string(argv[2]);
    if (argc>3) S3 = std::string(argv[3]);

    // The list of template geometries has been requested.

    if (Contains(S1,std::string("Templates")))
    {
        std::cout << "\n";
        std::cout << "The following template geometries are currently available for use in BEMUse:\n\n";
        std::cout << "Ellipsoid \n";
        std::cout << "SemiEllipsoid \n";
        std::cout << "Cylinder \n";
        std::cout << "Barge \n";
        std::cout << "TaperedSparBuoy \n";
        std::cout << "OC3SparBuoy - A reference spar geometry for FOWT platforms \n";
        std::cout << "TripleSpar - A three-legged spar type platform without a central column \n";
        std::cout << "OC4Semisub - A semisubmersible FOWT platform geometry  \n";
        std::cout << "WigleyHull - A standard reference ship hull geometry \n\n";
        std::cout << "If you would like your geometry on this list, please get in touch! \n";

        return true;
    }

    // Input information regarding a specific geometry has been requested

    if (Contains(S1,std::string("Inputs")))
    {
        // External files are being used to generate geometry....

        std::cout << "\n";

        if      (Contains(S2,std::string(".gdf")))
        {
            std::cout << "A WAMIT-style [.gdf] file has been specified as the input geometry.\n";
            std::cout << "Any [.bemin] dimension or discretisation input files will be ignored.\n";
        }
        else if (Contains(S2,std::string(".mar")))
        {
            std::cout << "A NEMOH-style [.mar] file has been specified as the input geometry.\n";
            std::cout << "Any [.bemin] dimension or discretisation input files will be ignored.\n";
        }
        else if (Contains(S2,std::string(".stl")))
        {
            std::cout << "An [.stl] file has been specified as the input geometry.\n";
            std::cout << "Any [.bemin] dimension or discretisation input files will be ignored.\n";
        }

        // Template geometries are being used...

        else if (Contains(S2,std::string("SemiEllipsoid")))
        {
            std::cout << "A submerged semi-ellipsoid geometry has been selected. If you used your imagination... this might be a hemisphere ;)\n\n";
            std::cout << "Three [BDRY] dimensional inputs are required. These specify the semi-axis of the ellipsoid.\n";
            std::cout << "One [EFS] dimensionalinput is required. This specifies the radius of the exterior free surface grid.\n\n";
            std::cout << "Two [BDRY] discretisation inputs are required. These specify the axial and azimuthal resolution of the grid.\n";
            std::cout << "One [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface grid.\n";
            std::cout << "Two [EFS] discretisation inputs are required. These specify the azimuthal and radial radial resolution of the exterior free surface grid.\n";
        }

        else if (Contains(S2,std::string("Ellipsoid")))
        {
            std::cout << "A submerged ellipsoid geometry has been selected.\n\n";
            std::cout << "Three [BDRY] dimensional inputs are required. These specify the semi-axis of the ellipsoid.\n\n";
            std::cout << "Two [BDRY] discretisation inputs are required. These specify the axial and azimuthal resolution of the grid.\n";
        }

        else if (Contains(S2,std::string("Cylinder")))
        {
            std::cout << "A submerged Half_Cylinder geometry has been selected.\n\n";
            std::cout << "Two [BDRY] dimensional inputs are required. These specify the radius and depth of the cylinder.\n";
            std::cout << "One [EFS] dimensional input is required. This specifies the radius of the exterior free surface grid.\n\n";
            std::cout << "Three [BDRY] discretisation inputs are required. These specify the azimuthal and radial resolution of the cylinder base and the vertical resolution of the cylinder.\n";
            std::cout << "One [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface.\n";
            std::cout << "One [EFS] discretisation inputs are required.This specifies the radial resolution of the exterior free surface.\n\n";
        }

        else if (Contains(S2,std::string("Barge")))
        {
            std::cout << "A Barge geometry has been selected.\n\n";
            std::cout << "Three [BDRY] dimensional inputs are required. These specify the length (X), width (Y) and draft (Z) of the barge.\n";
            std::cout << "Two [EFS] dimensional inputs are required. These specify the length (X) and width (Y) of the exterior free surface grid.\n\n";
            std::cout << "Three [BDRY] discretisation inputs are required. These specify the lengthwise, spanwise and vertical discretisation of the barge walls.\n";
    //        std::cout << "No [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface.\n";
            std::cout << "Two [EFS] discretisation inputs are required. These specify the lengthwise and spanwise resolution of the exterior free surface grid.\n";
        }

        else if (Contains(S2,std::string("TaperedSparBuoy")))
        {
            std::cout << "A Tapered Spar-Buoy geometry has been selected.\n\n";
            std::cout << "Five [BDRY] dimensional inputs are required. These specify the lower and upper radius of the cylinder sections, the draft of the buoy and the heights of the bottom and middle buoy sections.\n";
            std::cout << "One [EFS] dimensional input is required. This specifies the radius of the exterior free surface grid.\n\n";
            std::cout << "Five [BDRY] discretisation inputs are required. These specify:\n";
            std::cout << "i)   Azimuthal resolution.\n";
            std::cout << "ii)  Radial resolution of the base.\n";
            std::cout << "iii) Vertical resolution of the lower buoy section.\n";
            std::cout << "iv)  Vertical resolution of the middle buoy section.\n";
            std::cout << "v)   Vertical resolution of the upper buoy section.\n";
            std::cout << "One [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface grid.\n";
            std::cout << "One [EFS] discretisation input is required. This specifies the radial resolution of the exterior free surface grid.\n\n";
        }

        else if (Contains(S2,std::string("OC3SparBuoy")))
        {
            std::cout << "The OC3 floating offshore wind turbine platform geometry has been selected. Ref: doi:10.2172/979456 \n\n";
            std::cout << "One [BDRY] dimensional input is required. This parameter, 'S', scales the entire platform geometry.\n";
            std::cout << "A value of 1.0 returns the original OC3 geometry with the CoG_def and Mass_def as defined in the reference doc.\n";
            std::cout << "The linear position of the CoG is taken as S*COG_def.\n";
            std::cout << "The mass and translation inertial matrix parameters (1,1)- (3,3) are scaled by S^3.\n";
             std::cout << "The rotational inertial matrix parameters (4,4)- (6,6) are scaled by S^5.\n";
            std::cout << "One [EFS] dimensional input is required. This specifies the radius of the exterior free surface grid.\n\n";
            std::cout << "Five [BDRY] discretisation inputs are required. These specify:\n";
            std::cout << "i)   Azimuthal resolution.\n";
            std::cout << "ii)  Radial resolution of the base.\n";
            std::cout << "iii) Vertical resolution of the lower buoy section.\n";
            std::cout << "iv)  Vertical resolution of the middle buoy section.\n";
            std::cout << "v)   Vertical resolution of the upper buoy section.\n";
            std::cout << "One [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface grid.\n";
            std::cout << "One [EFS] discretisation input is required. This specifies the radial resolution of the exterior free surface grid.\n\n";
        }

        else if (Contains(S2,std::string("TripleSpar")))
        {
            std::cout << "A Hywind style triple-spar FOWT platform geometry has been selected.\n\n";
            std::cout << "Three [BDRY] dimensional inputs are required. These specify:\n";
            std::cout << "i)   Pontoon radius.\n";
            std::cout << "ii)  Pontoon depth.\n";
            std::cout << "iii) Radius of floater (refers to radius to the centre of each pontoon).\n";
    //        std::cout << "One [EFS] dimensional input is required. This specifies the radius of the exterior free surface grid.\n\n";
            std::cout << "Two [BDRY] discretisation inputs are required. These specify:\n";
            std::cout << "i)   Azimuthal resolution.\n";
            std::cout << "ii)  Vertical resolution of the pontoons.\n";
            std::cout << "One [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface grids.\n";
    //        std::cout << "One [EFS] discretisation input is required. This specifies the radial resolution of the exterior free surface grid.\n\n";
        }

        else if (Contains(S2,std::string("OC4Semisub")))
        {
            std::cout << "An OC4 style triple-spar FOWT platform geometry has been selected. Ref: doi:10.2172/1155123 \n\n";
            std::cout << "One [BDRY] dimensional input is required. This parameter, 'S', scales the entire platform geometry.\n";
            std::cout << "A value of 1.0 returns the original OC4 geometry with the CoG_def and Mass_def as defined in the reference doc.\n";
            std::cout << "The linear position of the CoG is taken as S*COG_def.\n";
            std::cout << "The mass and translation inertial matrix parameters (1,1)- (3,3) are scaled by S^3.\n";
             std::cout << "The rotational inertial matrix parameters (4,4)- (6,6) are scaled by S^5.\n\n";
    //        std::cout << "One [EFS] dimensional input is required. This specifies the radius of the exterior free surface grid.\n\n";
            std::cout << "Five [BDRY] discretisation inputs are required. These specify:\n";
            std::cout << "i)   Azimuthal resolution of pontoon.\n";
            std::cout << "ii)  Radial resolution of the pontoon base.\n";
            std::cout << "iii) Vertical resolution of the lower (heave plate) section of the pontoon.\n";
            std::cout << "iv)  Radial resolution of the upper (heave plate) section of the pontoon.\n";
            std::cout << "v)   Vertical resolution of the pontoon.\n";
            std::cout << "Nb: Parameters i), ii) and v) are used to define the resolution of the central column.\n";
            std::cout << "One [IFS] discretisation input is required. This specifies the radial resolution of the interior free surface grids.\n\n";
    //        std::cout << "One [EFS] discretisation input is required. This specifies the radial resolution of the exterior free surface grid.\n\n";
        }

        else if (Contains(S2,std::string("WigleyHull")))
        {
            std::cout << "An Wigley Hull geometry has been selected. \n\n";
            std::cout << "Three [BDRY] dimensional inputs are required. These specify the length (X), width (Y) and draft (Z) of the hull. \n\n";
            std::cout << "Two [BDRY] discretisation inputs are required. These specify the lengthwise and vertical resolution of the hull. \n\n";
        }

        else
        {
            std::cout << S2 << ": geometry type unknown. Please specify a valid geometry type.\n";

        }

        return true;

//        std::cout << "An X geometry has been selected.\n\n";
//        std::cout << "X [BDRY] dimensional inputs are required. These specify blah.\n";
//        std::cout << "X [EFS] dimensional input is required. These specify blah.\n\n";
//        std::cout << "X [BDRY] discretisation inputs are required. These specify blah.\n";
//        std::cout << "X [IFS] discretisation input is required. These specify blah.\n";
//        std::cout << "X [EFS] discretisation inputs are required.  These specify blah.\n\n";
    }

    //--- None of the enquiry commands have been fulfilled...
    //--- An analysis will be carried out.
    return false;
}

#endif // BEMUSE_INQUIRY_H
