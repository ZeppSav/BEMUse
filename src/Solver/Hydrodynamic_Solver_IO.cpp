//---------------------------------------------------------------
//------------Hydro solver Input/Output Functions----------------
//---------------------------------------------------------------

#include "../BEMUse_IO.h"
#include "Hydrodynamic_Solver.h"

namespace BEMUse
{

//--- Output functions

void Hydrodynamic_Radiation_Solver::Generate_Output_File(Boundary *B)
{
    // This creates the output files and fills preliminary data for the boundary / simulation.
    CreateDirectory(OutputDirectory);   // Create directory if it doesn't yet exist
    CreateDirectory(OutputPath);        // Create output file path if it doesn't yet exist
    Generate_Output_File_BEMUse(B);     // Generates output in BEMUse format
    Generate_Output_File_WAMIT(B);      // Generates output in WAMIT format
}

void Hydrodynamic_Radiation_Solver::Update_Output_File()
{
    // Cals both update functions
    Update_Output_File_BEMUse();
    Update_Output_File_WAMIT();
}

//--- Output (BEMUse Format)

void Hydrodynamic_Radiation_Solver::Generate_Output_File_BEMUse(Boundary *B)
{
    // The BEMUse output files is tailored to simplicity and succintness.
    // Furthermore, an attempt is made to store all data in a single file to simplify postprocessing.

    std::string FilePath = OutputPath + ".bem";
    std::ofstream file;

    file.open(FilePath, std::ofstream::out | std::ofstream::trunc); // Clear!
    file.close();

    // Write initial inputs
    file.open(FilePath);
    if (file.is_open())
    {
        file << "---------------------------------------------------------------\n";
        file << "------------------------ BEMUse Output ------------------------\n";
        file << "---------------------------------------------------------------\n";
        file << "\n";
        file << "------------ Environmental Parameters ------------\n";
        file << "Gravity:                 " << Gravity << " m s^{-2}.\n";
        file << "Density of Water:        " << Rho_wat << " kg m^{-3}.\n";
        file << "Atmospheric Pressure:    " << P_atm << " Pa.\n";
        if (Depth==Infinite) file << "Water Depth:             Infinite.\n";
        if (Depth==Finite) file   << "Water Depth:          "<< H << " m.\n";

        file << "\n";
        file << "------------ Wave Parameters ------------\n";
//        file << "Wave Amplitude:               "<< A << " m.\n";
//        file << "Incoming wave angles: (deg)        \n";
        std::string WAngles = "Incoming wave angles: (deg)        ";
//        WAngles.push_back("Incoming wave angles: (deg)        ");       // << std::endl;
        for (int i=0; i<BetaArray.size(); i++)   WAngles.append(PFS(BetaArray[i]*R2D));
        file << WAngles << std::endl;

        file << "\n";
        file << "------------ Grid Summary ------------\n";
//        std::string S = B->FileInfo;
//        if (S.empty())  file << "The geometry has been generated within BEMUse.\n";
        file << B->FileInfo << std::endl;
//        if (Geo_Vars->Geo==BEMUSE)  file << "The geometry is being generated within BEMUse based on the geometry defined in Geometry.dat.\n";
//        if (Geo_Vars->Geo==WAMIT)   file << "Geometry imported from .gdf file: " << IO_Vars->Input_Geo_File <<  std::endl;
//        if (Geo_Vars->Geo==NEMOH)   file << "Geometry imported from .mar/.dat file: " << IO_Vars->Input_Geo_File <<  std::endl;
        file << "Input Mesh has:          " << B->N_Nodes()       << " nodes.\n";
        file << "Input Mesh has:          " << B->N_Elements()    << " panels.\n";
        file << "Input Mesh has:          " << B->N_Aux_Elements()<< " auxiliary panels.\n";

        file << "\n";
        file << "------------ Geometry Summary ------------\n";
        file << "Nb: This is calculated by BEMUse with the approximation of flat surface panels."<< std::endl;
        file << "With increasing mesh refinement the values are expected to approach those of the continuous geometry.\n";
        file << "Volume:                  " << B->Volume << " m^3.\n";
        file << "Surface area:            " << B->Surface_Area << " m^2.\n";
        file << "Mass:                    " << B->Mass << " kg. (Assuming not pre-specified),\n";
        file << "Centre of Mass:          ["  << B->Centre_Gravity(0) << ", "
                                                << B->Centre_Gravity(1) << ", "
                                                << B->Centre_Gravity(2) << "] m.\n";
        file << "Centre of Buoyancy:      ["  << B->Centre_Buoyancy(0) << ", "
                                                << B->Centre_Buoyancy(1) << ", "
                                                << B->Centre_Buoyancy(2) << "] m.\n";

        file << "\n";
        file << "------------ Inertial Summary ------------\n";
//        file << "The mass matrix was calculated assuming uniform density within the geometry.\n";
//        file << "The entries M(4,4)-M(6,6) were calculated based on the Radii of Gyration specified.\n";
//        if (Inert_Vars->MassMatOption==2) file << "The mass matrix was imported from the file Inertia.dat."<< std::endl;
//        if (Inert_Vars->MassMatOption==1) file << "The mass matrix was calculated assuming uniform density within the geometry.\n";
//        if (Inert_Vars->MassMatOption==1) file << "The entries M(4,4)-M(6,6) were calculated based on the Radii of Gyration specified in the file Inertia.dat.\n";
//        if (Inert_Vars->MassMatOption==2) file << "The mass matrix was imported from the file Inertia.dat."<< std::endl;
        file << "\n";
        file  << "Mass Matrix: \n";
        file  << "|" <<  PFS(M_ij(0,0)) <<  PFS(M_ij(0,1)) <<  PFS(M_ij(0,2)) <<  PFS(M_ij(0,3)) <<  PFS(M_ij(0,4)) <<  PFS(M_ij(0,5)) << " |\n";
        file  << "|" <<  PFS(M_ij(1,0)) <<  PFS(M_ij(1,1)) <<  PFS(M_ij(1,2)) <<  PFS(M_ij(1,3)) <<  PFS(M_ij(1,4)) <<  PFS(M_ij(1,5)) << " |\n";
        file  << "|" <<  PFS(M_ij(2,0)) <<  PFS(M_ij(2,1)) <<  PFS(M_ij(2,2)) <<  PFS(M_ij(2,3)) <<  PFS(M_ij(2,4)) <<  PFS(M_ij(2,5)) << " |\n";
        file  << "|" <<  PFS(M_ij(3,0)) <<  PFS(M_ij(3,1)) <<  PFS(M_ij(3,2)) <<  PFS(M_ij(3,3)) <<  PFS(M_ij(3,4)) <<  PFS(M_ij(3,5)) << " |\n";
        file  << "|" <<  PFS(M_ij(4,0)) <<  PFS(M_ij(4,1)) <<  PFS(M_ij(4,2)) <<  PFS(M_ij(4,3)) <<  PFS(M_ij(4,4)) <<  PFS(M_ij(4,5)) << " |\n";
        file  << "|" <<  PFS(M_ij(5,0)) <<  PFS(M_ij(5,1)) <<  PFS(M_ij(5,2)) <<  PFS(M_ij(5,3)) <<  PFS(M_ij(5,4)) <<  PFS(M_ij(5,5)) << " |\n";
        file << "\n";
//        if (Inert_Vars->HydroMatOption==1)
//        {
//        file << "The hydrostatic stiffness matrix was calculated based on floater geometry and inertial parameters.\n";
//        file << "\n";
        file  << "Hydrostatic stiffness Matrix: \n";
        file  << "|" <<  PFS(C_ij(0,0)) <<  PFS(C_ij(0,1)) <<  PFS(C_ij(0,2)) <<  PFS(C_ij(0,3)) <<  PFS(C_ij(0,4)) <<  PFS(C_ij(0,5)) << " |\n";
        file  << "|" <<  PFS(C_ij(1,0)) <<  PFS(C_ij(1,1)) <<  PFS(C_ij(1,2)) <<  PFS(C_ij(1,3)) <<  PFS(C_ij(1,4)) <<  PFS(C_ij(1,5)) << " |\n";
        file  << "|" <<  PFS(C_ij(2,0)) <<  PFS(C_ij(2,1)) <<  PFS(C_ij(2,2)) <<  PFS(C_ij(2,3)) <<  PFS(C_ij(2,4)) <<  PFS(C_ij(2,5)) << " |\n";
        file  << "|" <<  PFS(C_ij(3,0)) <<  PFS(C_ij(3,1)) <<  PFS(C_ij(3,2)) <<  PFS(C_ij(3,3)) <<  PFS(C_ij(3,4)) <<  PFS(C_ij(3,5)) << " |\n";
        file  << "|" <<  PFS(C_ij(4,0)) <<  PFS(C_ij(4,1)) <<  PFS(C_ij(4,2)) <<  PFS(C_ij(4,3)) <<  PFS(C_ij(4,4)) <<  PFS(C_ij(4,5)) << " |\n";
        file  << "|" <<  PFS(C_ij(5,0)) <<  PFS(C_ij(5,1)) <<  PFS(C_ij(5,2)) <<  PFS(C_ij(5,3)) <<  PFS(C_ij(5,4)) <<  PFS(C_ij(5,5)) << " |\n";
        file << "\n";
//        file  << "Hydrostatic stiffness Matrix (Without gravity loads): \n";
//        file  << "|" <<  PFS(C_ij_NM(0,0)) <<  PFS(C_ij_NM(0,1)) <<  PFS(C_ij_NM(0,2)) <<  PFS(C_ij_NM(0,3)) <<  PFS(C_ij_NM(0,4)) <<  PFS(C_ij_NM(0,5)) << " |\n";
//        file  << "|" <<  PFS(C_ij_NM(1,0)) <<  PFS(C_ij_NM(1,1)) <<  PFS(C_ij_NM(1,2)) <<  PFS(C_ij_NM(1,3)) <<  PFS(C_ij_NM(1,4)) <<  PFS(C_ij_NM(1,5)) << " |\n";
//        file  << "|" <<  PFS(C_ij_NM(2,0)) <<  PFS(C_ij_NM(2,1)) <<  PFS(C_ij_NM(2,2)) <<  PFS(C_ij_NM(2,3)) <<  PFS(C_ij_NM(2,4)) <<  PFS(C_ij_NM(2,5)) << " |\n";
//        file  << "|" <<  PFS(C_ij_NM(3,0)) <<  PFS(C_ij_NM(3,1)) <<  PFS(C_ij_NM(3,2)) <<  PFS(C_ij_NM(3,3)) <<  PFS(C_ij_NM(3,4)) <<  PFS(C_ij_NM(3,5)) << " |\n";
//        file  << "|" <<  PFS(C_ij_NM(4,0)) <<  PFS(C_ij_NM(4,1)) <<  PFS(C_ij_NM(4,2)) <<  PFS(C_ij_NM(4,3)) <<  PFS(C_ij_NM(4,4)) <<  PFS(C_ij_NM(4,5)) << " |\n";
//        file  << "|" <<  PFS(C_ij_NM(5,0)) <<  PFS(C_ij_NM(5,1)) <<  PFS(C_ij_NM(5,2)) <<  PFS(C_ij_NM(5,3)) <<  PFS(C_ij_NM(5,4)) <<  PFS(C_ij_NM(5,5)) << " |\n";
//        file << "\n";
//        }
//        if (Inert_Vars->HydroMatOption==2)
//        {
//            file << "The hydrostatic stiffness matrix was imported from the file Inertia.dat."<< std::endl;

//            file  << "Hydrostatic stiffness Matrix: \n";
//            file  << "|" <<  PFS(C_ij(0,0)) <<  PFS(C_ij(0,1)) <<  PFS(C_ij(0,2)) <<  PFS(C_ij(0,3)) <<  PFS(C_ij(0,4)) <<  PFS(C_ij(0,5)) << " |\n";
//            file  << "|" <<  PFS(C_ij(1,0)) <<  PFS(C_ij(1,1)) <<  PFS(C_ij(1,2)) <<  PFS(C_ij(1,3)) <<  PFS(C_ij(1,4)) <<  PFS(C_ij(1,5)) << " |\n";
//            file  << "|" <<  PFS(C_ij(2,0)) <<  PFS(C_ij(2,1)) <<  PFS(C_ij(2,2)) <<  PFS(C_ij(2,3)) <<  PFS(C_ij(2,4)) <<  PFS(C_ij(2,5)) << " |\n";
//            file  << "|" <<  PFS(C_ij(3,0)) <<  PFS(C_ij(3,1)) <<  PFS(C_ij(3,2)) <<  PFS(C_ij(3,3)) <<  PFS(C_ij(3,4)) <<  PFS(C_ij(3,5)) << " |\n";
//            file  << "|" <<  PFS(C_ij(4,0)) <<  PFS(C_ij(4,1)) <<  PFS(C_ij(4,2)) <<  PFS(C_ij(4,3)) <<  PFS(C_ij(4,4)) <<  PFS(C_ij(4,5)) << " |\n";
//            file  << "|" <<  PFS(C_ij(5,0)) <<  PFS(C_ij(5,1)) <<  PFS(C_ij(5,2)) <<  PFS(C_ij(5,3)) <<  PFS(C_ij(5,4)) <<  PFS(C_ij(5,5)) << " |\n";
//            file << "\n";
//            file  << "Hydrostatic stiffness Matrix (Without gravity loads): Not calculated. Matrix has been imported (see above).\n";
//            file << "\n";
//        }
        file << "------------ Output of Hydrodynamic Analysis: ------------\n";
        file << "\n";
        file << "\n";
    }
}

void Hydrodynamic_Radiation_Solver::Update_Output_File_BEMUse()
{
    // This updates the output file using the BEMUse format
    std::string FilePath = OutputPath + ".bem";
    std::ofstream file;

    file.open(FilePath, std::ios_base::app);
    if (file.is_open())
    {
//        QTextStream stream(&file);
        file << "----------\n";
        file << "Frequency:       " << Omega << " s^{-1}."    << std::endl;
        file << "Period:          " << Period << " s."    << std::endl;
        file << "Wavenumber:      " << Kappa << " m^{-1}."    << std::endl;
        file << "\n";

//        std::string Title1;
//        D.append("Motion");
//        Title1.append()
//        PST()
        file << PST("Motion")   << PST(" A_x")      << PST(" B_x")      << PST(" A_y")      << PST(" B_y")      << PST(" A_z")      << PST(" B_z")
                                << PST(" A_pitch")  << PST(" B_pitch")  << PST(" A_roll")   << PST(" B_roll")   << PST(" A_yaw")    << PST(" B_yaw") << std::endl;
//        file << "Motion             A_x           B_x           A_y           B_y           A_z           B_z           A_pitch       B_pitch       A_roll        B_roll        A_yaw         B_yaw        "<<  std::endl;
        file << PST("Surge")    <<  PFS(A_ij(0,0)) << PFS(B_ij(0,0))  << PFS(A_ij(0,1)) << PFS(B_ij(0,1))  <<  PFS(A_ij(0,2)) << PFS(B_ij(0,2))
                                <<  PFS(A_ij(0,3)) << PFS(B_ij(0,3)) <<  PFS(A_ij(0,4)) << PFS(B_ij(0,4))  <<  PFS(A_ij(0,5)) << PFS(B_ij(0,5)) << std::endl;
//                                        i = ID_Sway;
        file << PST("Sway")     <<  PFS(A_ij(1,0)) << PFS(B_ij(1,0))  <<  PFS(A_ij(1,1)) << PFS(B_ij(1,1))  <<  PFS(A_ij(1,2)) << PFS(B_ij(1,2))
                                <<  PFS(A_ij(1,3)) << PFS(B_ij(1,3)) <<  PFS(A_ij(1,4)) << PFS(B_ij(1,4))  <<  PFS(A_ij(1,5)) << PFS(B_ij(1,5)) << std::endl;
//                                        i = ID_Heave;
        file << PST("Heave")    <<  PFS(A_ij(2,0)) << PFS(B_ij(2,0))  <<  PFS(A_ij(2,1)) << PFS(B_ij(2,1))  <<  PFS(A_ij(2,2)) << PFS(B_ij(2,2))
                                <<  PFS(A_ij(2,3)) << PFS(B_ij(2,3))  <<  PFS(A_ij(2,4)) << PFS(B_ij(2,4))  <<  PFS(A_ij(2,5)) << PFS(B_ij(2,5)) << std::endl;
//                                        i = ID_Roll;
        file << PST("Roll")     << PFS(A_ij(3,0)) << PFS(B_ij(3,0))  <<  PFS(A_ij(3,1)) << PFS(B_ij(3,1))  <<  PFS(A_ij(3,2)) << PFS(B_ij(3,2))
                                <<  PFS(A_ij(3,3)) << PFS(B_ij(3,3))  <<  PFS(A_ij(3,4)) << PFS(B_ij(3,4))  <<  PFS(A_ij(3,5)) << PFS(B_ij(3,5)) << std::endl;
//Matrix                                        i = ID_Pitch;
        file << PST("Pitch")    <<  PFS(A_ij(4,0)) << PFS(B_ij(4,0))  <<  PFS(A_ij(4,1)) << PFS(B_ij(4,1))  <<  PFS(A_ij(4,2)) << PFS(B_ij(4,2))
                                <<  PFS(A_ij(4,3)) << PFS(B_ij(4,3))  <<  PFS(A_ij(4,4)) << PFS(B_ij(4,4))  <<  PFS(A_ij(4,5)) << PFS(B_ij(4,5)) << std::endl;
//                                        i = ID_Yaw;
        file << PST("Yaw")      <<  PFS(A_ij(5,0)) << PFS(B_ij(5,0))  <<  PFS(A_ij(5,1)) << PFS(B_ij(5,1))  <<  PFS(A_ij(5,2)) << PFS(B_ij(5,2))
                                <<  PFS(A_ij(5,3)) << PFS(B_ij(5,3))  <<  PFS(A_ij(5,4)) << PFS(B_ij(5,4))  <<  PFS(A_ij(5,5)) << PFS(B_ij(5,5)) << std::endl;

        file << "\n";
        file << "\n";

        file << "Froude-Krylov Forces:\n";
        std::string WAs, HDF, FKX, FKY, FKZ, FKRX, FKRY, FKRZ;
        WAs.append(PST("Wave Heading: (deg) ",30));
        HDF.append(PST("Motion: ",30));
        FKX.append(PST("Surge: ",30));
        FKY.append(PST("Sway: ",30));
        FKZ.append(PST("Heave: ",30));
        FKRX.append(PST("Pitch: ",30));
        FKRY.append(PST("Roll: ",30));
        FKRZ.append(PST("Yaw: ",30));
        for (int i=0; i<BetaArray.size(); i++)
        {
            WAs.append(PFS(BetaArray[i]*R2D,30));
            HDF.append(PST(" Mag(Fh)",15));
            HDF.append(PST(" Arg(Fh)",15));
            FKX.append(PFS(std::abs(FK_i(0,i))/Gravity,15));
            FKX.append(PFS(std::arg(FK_i(0,i)),15));
            FKY.append(PFS(std::abs(FK_i(1,i))/Gravity,15));
            FKY.append(PFS(std::arg(FK_i(1,i)),15));
            FKZ.append(PFS(std::abs(FK_i(2,i))/Gravity,15));
            FKZ.append(PFS(std::arg(FK_i(2,i)),15));
            FKRX.append(PFS(std::abs(FK_i(3,i))/Gravity,15));
            FKRX.append(PFS(std::arg(FK_i(3,i)),15));
            FKRY.append(PFS(std::abs(FK_i(4,i))/Gravity,15));
            FKRY.append(PFS(std::arg(FK_i(4,i)),15));
            FKRZ.append(PFS(std::abs(FK_i(5,i))/Gravity,15));
            FKRZ.append(PFS(std::arg(FK_i(5,i)),15));
        }
        file << WAs << std::endl;
        file << HDF << std::endl;
        file << FKX << std::endl;
        file << FKY << std::endl;
        file << FKZ << std::endl;
        file << FKRX << std::endl;
        file << FKRY << std::endl;
        file << FKRZ << std::endl;

        file << "\n";
        file << "\n";

        file << "Excitation Forces:\n";
        std::string  HDE, EFX, EFY, EFZ, EFRX, EFRY, EFRZ;
//        QFillAppend(WAs,"Wave Heading: (deg) ",30);
        HDE.append(PST("Motion: ",30));
        EFX.append(PST("Surge: ",30));
        EFY.append(PST("Sway: ",30));
        EFZ.append(PST("Heave: ",30));
        EFRX.append(PST("Pitch: ",30));
        EFRY.append(PST("Roll: ",30));
        EFRZ.append(PST("Yaw: ",30));
        for (int i=0; i<BetaArray.size(); i++)
        {
//            WAs,Sim_Vars->Wave_Angles[i]*R2D,30));
            HDE.append(PST(" Mag(Xh)",15));
            HDE.append(PST(" Arg(Xh)",15));
            EFX.append(PFS(std::abs(F_i(0,i))/Gravity,15));
            EFX.append(PFS(std::arg(F_i(0,i)),15));
            EFY.append(PFS(std::abs(F_i(1,i))/Gravity,15));
            EFY.append(PFS(std::arg(F_i(1,i)),15));
            EFZ.append(PFS(std::abs(F_i(2,i))/Gravity,15));
            EFZ.append(PFS(std::arg(F_i(2,i)),15));
            EFRX.append(PFS(std::abs(F_i(3,i))/Gravity,15));
            EFRX.append(PFS(std::arg(F_i(3,i)),15));
            EFRY.append(PFS(std::abs(F_i(4,i))/Gravity,15));
            EFRY.append(PFS(std::arg(F_i(4,i)),15));
            EFRZ.append(PFS(std::abs(F_i(5,i))/Gravity,15));
            EFRZ.append(PFS(std::arg(F_i(5,i)),15));
        }
        file << WAs << std::endl;
        file << HDE << std::endl;
        file << EFX << std::endl;
        file << EFY << std::endl;
        file << EFZ << std::endl;
        file << EFRX << std::endl;
        file << EFRY << std::endl;
        file << EFRZ << std::endl;

        file << "\n";
        file << "\n";

//        file << "Response Amplitude Operator:\n";
//        std::string HDRAO, RAOX, RAOY, RAOZ, RAORX, RAORY, RAORZ;
////        QFillAppend(WAs,"Wave Heading: (deg) ",30);
//        QFillAppend(HDRAO,"Motion: ",30);
//        QFillAppend(RAOX,"Surge: ",30);
//        QFillAppend(RAOY,"Sway: ",30);
//        QFillAppend(RAOZ,"Heave: ",30);
//        QFillAppend(RAORX,"Pitch: ",30);
//        QFillAppend(RAORY,"Roll: ",30);
//        QFillAppend(RAORZ,"Yaw: ",30);
//        for (int i=0; i<BetaArray.size(); i++)
//        {
////            QFillAppend(WAs,Sim_Vars->Wave_Angles[i]*R2D,30);
//            QFillAppend(HDRAO,"Mag(Xh)",15);
//            QFillAppend(HDRAO,"Arg(Xh)",15);
//            QFillAppendLong(RAOX,std::abs(RAO_i(0,i)),15);
//            QFillAppendLong(RAOX,std::arg(RAO_i(0,i)),15);
//            QFillAppendLong(RAOY,std::abs(RAO_i(1,i)),15);
//            QFillAppendLong(RAOY,std::arg(RAO_i(1,i)),15);
//            QFillAppendLong(RAOZ,std::abs(RAO_i(2,i)),15);
//            QFillAppendLong(RAOZ,std::arg(RAO_i(2,i)),15);
//            QFillAppendLong(RAORX,std::abs(RAO_i(3,i)),15);
//            QFillAppendLong(RAORX,std::arg(RAO_i(3,i)),15);
//            QFillAppendLong(RAORY,std::abs(RAO_i(4,i)),15);
//            QFillAppendLong(RAORY,std::arg(RAO_i(4,i)),15);
//            QFillAppendLong(RAORZ,std::abs(RAO_i(5,i)),15);
//            QFillAppendLong(RAORZ,std::arg(RAO_i(5,i)),15);
//        }
//        file << WAs << std::endl;
//        file << HDRAO << std::endl;
//        file << RAOX << std::endl;
//        file << RAOY << std::endl;
//        file << RAOZ << std::endl;
//        file << RAORX << std::endl;
//        file << RAORY << std::endl;
//        file << RAORZ << std::endl;

        file << "----------" ;
        file << "\n";
        file << "\n";

        // Print sim time if completed

//        if (Sim_Vars->Timestep==Sim_Vars->NF-1)
//        {
//            file << "----------\n";
//            file << "Simulation summary:\n";
//            file << "----------\n";
//            Real Seconds = Real(Sim_Vars->Sim_Clock.restart())/1000.0;
//            file << "Simulation time: " << Seconds << " seconds.\n";
//        }

//        std::cout << Omega      <<" "<< A_ij(0,0) <<" "<< B_ij(0,0)
//                                <<" "<< A_ij(1,1) <<" "<< B_ij(1,1)
//                                <<" "<< A_ij(2,2) <<" "<< B_ij(2,2)
//                                <<" "<< A_ij(3,3) <<" "<< B_ij(3,3)
//                                <<" "<< A_ij(4,4) <<" "<< B_ij(4,4)
//                                <<" "<< A_ij(5,5) <<" "<< B_ij(5,5) << std::endl;

//        std::cout << Omega  <<" "<< std::abs(F_i(0)) <<" "<< std::arg(F_i(0))*R2D
//                            <<" "<< std::abs(F_i(1)) <<" "<< std::arg(F_i(1))*R2D
//                            <<" "<< std::abs(F_i(2)) <<" "<< std::arg(F_i(2))*R2D
//                            <<" "<< std::abs(F_i(3)) <<" "<< std::arg(F_i(3))*R2D
//                            <<" "<< std::abs(F_i(4)) <<" "<< std::arg(F_i(4))*R2D
//                            <<" "<< std::abs(F_i(5)) <<" "<< std::arg(F_i(5))*R2D << std::endl;

//        std::cout << Period <<" "<< std::abs(RAO_i(0)) <<" "<< std::arg(RAO_i(0))*R2D
//                            <<" "<< std::abs(RAO_i(1)) <<" "<< std::arg(RAO_i(1))*R2D
//                            <<" "<< std::abs(RAO_i(2)) <<" "<< std::arg(RAO_i(2))*R2D
//                            <<" "<< std::abs(RAO_i(3)) <<" "<< std::arg(RAO_i(3))*R2D
//                            <<" "<< std::abs(RAO_i(4)) <<" "<< std::arg(RAO_i(4))*R2D
//                            <<" "<< std::abs(RAO_i(5)) <<" "<< std::arg(RAO_i(5))*R2D << std::endl;


    }
}

//--- Output (WAMIT Format)

void Hydrodynamic_Radiation_Solver::Generate_Output_File_WAMIT(Boundary *B)
{
    // The WAMIT files are created and prepared. This is done in such a way to avoid any issues
    // A nice feature is prepared below:

    std::vector<std::string> Paths;
    Paths.push_back(OutputPath + ".1");     // Radiation force file
    Paths.push_back(OutputPath + ".2");     // Haskind force file
    Paths.push_back(OutputPath + ".2fk");   // Haskind force file (Froude-Krylov)
    Paths.push_back(OutputPath + ".2sc");   // Haskind force file (Scattering)
//    Paths.push_back(OutputPath + ".3");     // Diffraction force file
//    Paths.push_back(OutputPath + ".3fk");   // Diffraction force file (Froude-Krylov)
//    Paths.push_back(OutputPath + ".3sc");   // Diffraction force file (Scattering)
    Paths.push_back(OutputPath + ".4");     // Response amplitude operator (Haskind force)
    Paths.push_back(OutputPath + ".6");     // Wave elevations
    Paths.push_back(OutputPath + ".hst");   // Hydrostatic stiffness matrix
    Paths.push_back(OutputPath + ".mmx");   // Mass matrix
    Paths.push_back(OutputPath + "_KochinRad.bem");     // Kochin function radiation terms
    Paths.push_back(OutputPath + "_KochinDiff.bem");    // Kochin function diffraction terms

    std::ofstream file;

    // Clear files
    for (int i=0; i<Paths.size(); i++){
        file.open(Paths[i], std::ofstream::out | std::ofstream::trunc); // Clear!
        file.close();
    }

    // Store the hydrostatic stiffness matrix
    std::string HSTPath = OutputPath + ".hst";
    file.open(HSTPath);
    if (file.is_open())
    {
        for (int i=0; i<NDOF; i++){
            for (int j=0; j<NDOF; j++)  file << "     " << i+1 << "     " << j+1 << "   " << PFS(C_ij(i,j)/Rho/Gravity) << std::endl;
        }
    }
    file.close();

    // Store the mass matrix
    std::string MMXPath = OutputPath + ".mmx";
    file.open(MMXPath);
    if (file.is_open())
    {
        file << std::endl;
        file << "Gravity:     "  << Grav << "     Length scale:        1.00000"  << std::endl;
        file << std::endl;
        file << "Volume     "  << B->Volume << std::endl;
        file << "Center of Buoyancy (Xb,Yb,Zb):   "  << PFS(B->Centre_Buoyancy(0))  << PFS(B->Centre_Buoyancy(1))   << PFS(B->Centre_Buoyancy(2)) << std::endl;
        file << "Center of Gravity  (Xg,Yg,Zg):   "  << PFS(B->Centre_Gravity(0))   << PFS(B->Centre_Gravity(1))    << PFS(B->Centre_Gravity(2)) << std::endl;
        file << std::endl;
        for (int i=0; i<NDOF; i++){
            for (int j=0; j<NDOF; j++)  file << "     " << i+1 << "     " << j+1 << "   " << PFS(M_ij(i,j)/Rho) << std::endl;
        }
    }
    file.close();

    // Specify free surface wave elevation output files
    if (!Wave_Nodes.empty())    // No external nodes, ignore!
    {
        std::vector<Vector3> WaveNodePos;
        for (SP_Node N: Wave_Nodes) WaveNodePos.push_back(N->Position_Global());
        std::string WEPath = OutputPath + "_WaveElPos.bem";
        file.open(WEPath);
        if (file.is_open()){
            for (int i=0; i<WaveNodePos.size(); i++) file << i << "     " << PFS(WaveNodePos[i](0)) << PFS(WaveNodePos[i](1)) << std::endl;
        }

    }

}

void Hydrodynamic_Radiation_Solver::Update_Output_File_WAMIT()
{
    // This function updates the output files in WAMIT style.
    std::ofstream file;

    // Update the radiation forces
    std::string RadPath = OutputPath + ".1";
    file.open(RadPath, std::ios_base::app);
    if (file.is_open())
    {
     for (int i=0; i<NDOF; i++){
         for (int j=0; j<NDOF; j++){
             if      (Omega==0.0)           file << PFS(0.0)        << i+1 << "     " << j+1 << "     " << PFS(std::real(A_ij(i,j))) << std::endl;
             else if (Frequency==FreqInf)   file << PFS(-1.0)       << i+1 << "     " << j+1 << "     " << PFS(std::real(A_ij(i,j))) << std::endl;
             else                           file << PFS(Frequency)  << i+1 << "     " << j+1 << "     " << PFS(A_ij(i,j)) << PFS(B_ij(i,j)) << std::endl;
         }
     }
    }
    file.close();

    // Update the excitation forces
    std::string DiffPath = OutputPath + ".2";
    file.open(DiffPath, std::ios_base::app);
    if (file.is_open())
    {
        for (int j=0; j<NBeta; j++){
            for (int i=0; i<NDOF; i++){
                if (Omega==0.0)             continue;       // No wave!
                if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                file << PFS(Frequency)  << PFS(BetaArray[j]*R2D) << i+1 << "  "
                                        << PFS(std::abs(F_i(i,j))/Gravity)  << PFS(std::arg(F_i(i,j)))
                                        << PFS(std::real(F_i(i,j))/Gravity) << PFS(std::imag(F_i(i,j))/Gravity)<< std::endl;
            }
        }
    }
    file.close();

    // Update the Haskind Froude-Krylov forces
    std::string DiffPathFK = OutputPath + ".2fk";
    file.open(DiffPathFK,std::ios_base::app);
    if (file.is_open())
    {
        for (int j=0; j<NBeta; j++){
            for (int i=0; i<NDOF; i++){
                if (Omega==0.0)             continue;       // No wave!
                if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                file << PFS(Frequency) << PFS(BetaArray[j]*R2D) << i+1 << "  "
                                    << PFS(std::abs(FK_i(i,j))/Gravity) <<  PFS(std::arg(FK_i(i,j)))
                                    << PFS(std::real(FK_i(i,j))/Gravity)<< PFS(std::imag(FK_i(i,j))/Gravity) << std::endl;
            }
        }
    }
    file.close();

    // Update the Haskind scattering forces
    std::string DiffPathSC = OutputPath + ".2sc";
    file.open(DiffPathSC,std::ios_base::app);
    if (file.is_open())
    {
        for (int j=0; j<NBeta; j++){
            for (int i=0; i<NDOF; i++){
                if (Omega==0.0)             continue;       // No wave!
                if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                file << PFS(Frequency)    << PFS(BetaArray[j]*R2D) << i+1 << "  "
                                         << PFS(std::abs(SC_i(i,j))/Gravity) <<  PFS(-std::arg(SC_i(i,j)))
                                         << PFS(std::real(SC_i(i,j))/Gravity) << PFS(-std::imag(SC_i(i,j))/Gravity)<< std::endl;
            }
        }
    }
    file.close();

    // Update the response amplitude operator terms
    std::string DiffPathRAO = OutputPath + ".4";
    file.open(DiffPathRAO,std::ios_base::app);
    if (file.is_open())
    {
        for (int j=0; j<NBeta; j++){
            for (int i=0; i<NDOF; i++){
                if (Omega==0.0)             continue;       // No wave!
                if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                file << PFS(Frequency) << PFS(BetaArray[j]*R2D) << i+1 << "  "   << PFS(std::abs(RAO_i(i,j))) << PFS(std::arg(RAO_i(i,j)))
                                                                             << PFS(std::real(RAO_i(i,j))) << PFS(std::imag(RAO_i(i,j)))<< std::endl;
            }
        }
    }
    file.close();

    // Update the kochin function for the radiation potential
    std::string DiffPathKochRad = OutputPath + "_KochinRad.bem";
    file.open(DiffPathKochRad,std::ios_base::app);
    if (file.is_open())
    {
        for (int i=0; i<NKoch; i++){
            Real Theta = (0.5+1.0*i)*TwoPI/NKoch;
            for (int d=0; d<NDOF; d++){
                if (Omega==0.0)             continue;       // No wave!
                if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                file    << PFS(Frequency) << d+1 << "  " << PFS(Theta*R2D)
                        << PFS(std::abs(KochinRad(i,d)))     << PFS(std::arg(KochinRad(i,d)))
                        << PFS(std::real(KochinRad(i,d)))   << PFS(std::imag(KochinRad(i,d))) << std::endl;
            }
        }
    }
    file.close();

    // Update the kochin function for the radiation potential
    std::string DiffPathKochDiff = OutputPath + "_KochinDiff.bem";
    file.open(DiffPathKochDiff,std::ios_base::app);
    if (file.is_open())
    {
        for (int d=0; d<NBeta; d++){
            for (int i=0; i<NKoch; i++){
                 Real Theta = (0.5+1.0*i)*TwoPI/NKoch;
                 if (Omega==0.0)             continue;       // No wave!
                 if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                 file  << PFS(Frequency) << PFS(BetaArray[d]*R2D) << "  " << PFS(Theta*R2D)
                         << PFS(std::abs(KochinDiff(i,d)))     << PFS(std::arg(KochinDiff(i,d)))
                         << PFS(std::real(KochinDiff(i,d)))   << PFS(std::imag(KochinDiff(i,d))) << std::endl;
            }
        }
    }
    file.close();

    // Specify free surface wave elevation output files
    if (!Wave_Nodes.empty())    // No external nodes, ignore!
    {
        std::string WERadPath = OutputPath + ".6r";
        file.open(WERadPath,std::ios_base::app);
        if (file.is_open())
        {
            for (int i=0; i<Wave_Nodes.size(); i++){
                if (Omega==0.0)             continue;       // No wave!
                if (Frequency==FreqInf)     continue;       // Crazy wave train (Blizzard of Ozz!!)
                file   << PFS(Frequency) << i << "     "
                     << PFS(std::abs(ExtRadMat(i,0)))     << PFS(std::arg(ExtRadMat(i,0)))
                     << PFS(std::abs(ExtRadMat(i,1)))     << PFS(std::arg(ExtRadMat(i,1)))
                     << PFS(std::abs(ExtRadMat(i,2)))     << PFS(std::arg(ExtRadMat(i,2)))
                     << PFS(std::abs(ExtRadMat(i,3)))     << PFS(std::arg(ExtRadMat(i,3)))
                     << PFS(std::abs(ExtRadMat(i,4)))     << PFS(std::arg(ExtRadMat(i,4)))
                     << PFS(std::abs(ExtRadMat(i,5)))     << PFS(std::arg(ExtRadMat(i,5))) << std::endl;
            }
        }
    }
}

}
