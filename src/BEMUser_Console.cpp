//---------------------------------------------------------------
//------------------ BEMuse Console Functions--------------------
//---------------------------------------------------------------

#include "BEMUse_IO.h"
#include "BEMUser_Console.h"

//--- Geometry types
#include "Boundary/Ellipsoid.h"
#include "Boundary/Volume_of_Revolution.h"
#include "Boundary/Triple_Spar.h"
#include "Boundary/Ship_Hulls.h"
#include "Boundary/STL_Geo.h"
#include "Boundary/GDF_Geo.h"
#include "Boundary/MAR_Geo.h"
#include "Boundary/FOWT_Platforms.h"
#include "Boundary/Barge.h"

//--- Solver params
// #include "omp.h"

//--- Solver types
#include "Solver/Aerodynamic_Solver.h"
#include "Solver/Hydrodynamic_Solver.h"

//--- Geometry spec

void BEMUse_Console::Specify_Geometry(int argc, char *argv[])
{
    if (argc<2) {
        std::cout << "Geometry type has not been specified.\n";
        return;
    }

    // Specify geometry type
    bool GeoImported = false;
    std::string S = std::string(argv[1]);

    // Imported geometries
    if      (Contains(S,std::string(".gdf")))               {Boundary = new BEMUse::GDF_Geometry(); GeoImported = true;}
    else if (Contains(S,std::string(".mar")))               {Boundary = new BEMUse::MAR_Geometry(); GeoImported = true;}
    else if (Contains(S,std::string(".stl")))               {Boundary = new BEMUse::STL_Geometry(); GeoImported = true;}

    // Template geometries
    else if (Contains(S,std::string("SemiEllipsoid")))      Boundary = new BEMUse::Semi_Ellipsoid();
    else if (Contains(S,std::string("Ellipsoid")))          Boundary = new BEMUse::Ellipsoid();
    else if (Contains(S,std::string("Cylinder")))           Boundary = new BEMUse::Half_Cylinder();
    else if (Contains(S,std::string("Barge")))              Boundary = new BEMUse::Barge();
    else if (Contains(S,std::string("TaperedSparBuoy")))    Boundary = new BEMUse::Tapered_SparBuoy();
    else if (Contains(S,std::string("OC3SparBuoy")))        Boundary = new BEMUse::OC3_SparBuoy();
    else if (Contains(S,std::string("TripleSpar")))         Boundary = new BEMUse::Triple_Spar();
    else if (Contains(S,std::string("OC4Semisub")))         Boundary = new BEMUse::SemiSub_OC4();
    else if (Contains(S,std::string("WigleyHull")))         Boundary = new BEMUse::Wigley_Hull();
    else
    {
        std::cout << S << ": geometry type unknown. Please specify a valid geometry type.\n";
        return;
    }

    //--- Setup of imported geometry
    if  (GeoImported){
        bool ImportSuccess = Boundary->Read_Input_File(S);
        if (ImportSuccess)  Boundary->Setup();
        else{
            delete Boundary;
            Boundary = NULL;
        }
        return;
    }

    //--- ifstream vars

    std::string FilePath;
    std::string line;

    //--- Setup of template geometry

    // Specify geometry dimensions. These are read in from a file named "Dimensions.bemin"
    std::vector<Real> Dim, DimAux, DimExt;
    FilePath = std::string("Input") + '/' + std::string("Dimensions.bemin");
    std::ifstream Dimfile;
    Dimfile.open(FilePath);
    if (Dimfile.is_open())
    {
        while ( std::getline(Dimfile,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields[0] == "BDRY")    Dim.push_back(std::stod(Fields[1]));
            if (Fields[0] == "IFS")     DimAux.push_back(std::stod(Fields[1]));
            if (Fields[0] == "EFS")     DimExt.push_back(std::stod(Fields[1]));
        }
    }
    Dimfile.close();

    // Specify geometry discretisation. These are read in from a file named "Dimensions.bemin"
    std::vector<int> Disc, DiscAux, DiscExt;
    FilePath = std::string("Input") + '/' + std::string("Discretisation.bemin");
    std::ifstream Discfile;
    Discfile.open(FilePath);
    if (Discfile.is_open())
    {
        while ( std::getline(Discfile,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields[0] == "BDRY")    Disc.push_back(std::stoi(Fields[1]));
            if (Fields[0] == "IFS")     DiscAux.push_back(std::stoi(Fields[1]));
            if (Fields[0] == "EFS")     DiscExt.push_back(std::stoi(Fields[1]));
        }
    }
    Discfile.close();

    // Specify geometry discretisation. These are read in from a file named "Dimensions.bemin"
    std::vector<bool> Flags;
    FilePath = std::string("Input") + '/' + std::string("Flags.bemin");
    std::ifstream Flagfile;
    Flagfile.open(FilePath);
    if (Flagfile.is_open())
    {
        while ( std::getline(Flagfile,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields[0] == "TRUE")  Flags.push_back(true);
            if (Fields[0] == "FALSE") Flags.push_back(false);
        }
    }
    Flagfile.close();

    //--- Generate geometry
    Boundary->Setup();
}

//--- Solver parameters

void BEMUse_Console::Specify_Compute_Params(int argc, char *argv[])
{
    // Check computational parameters

    // Checks for specification of threads.
    int NTM = omp_get_max_threads();
    int NT = NTM;
    for (int i=1; i<argc; i++)
    {
        std::string S = std::string(argv[i]);
        if (Contains(S,std::string("-t"))){
            std::vector<std::string> Fields = Split(S,'t');     // Split line into segments
            NT = std::stoi(Fields[1]);
        }
    }
    omp_set_num_threads(NT);
    std::cout << "\nBEMUse has been specified to run on " << NT << " of " << NTM << " available threads.\n\n";
}

//--- Solver spec

void BEMUse_Console::Specify_Solver(int argc, char *argv[])
{
    // Collects solver information and then generates the solver object
    if (argc<3) {
        std::cout << "Output data names not specified.\n";
        return;
    }

    if (Boundary==NULL){
        std::cout << "Geometry not yet specified. Solver initialisation & preprocessing aborting.\n";
        return;
    }

    //--- ifstream vars

    std::string FilePath;
    std::string line;
    std::vector<BEMUse::Parameter> Params;

    // Specify frequency list. These are read in from a file named "Frequencies.bemin"
    FilePath = std::string("Input") + '/' + std::string("Frequencies.bemin");
    std::ifstream Freqfile;
    Freqfile.open(FilePath);
    if (Freqfile.is_open())
    {
        while ( std::getline(Freqfile,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            Frequencies.push_back(std::stod(Fields[0]));
        }
    }
    Freqfile.close();

    FilePath = std::string("Input") + '/' + std::string("WaveAngles.bemin");
    std::ifstream WAfile;
    WAfile.open(FilePath);
    if (WAfile.is_open())
    {
        while ( std::getline(WAfile,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            Params.push_back(BEMUse::Parameter("WaveAngle",Real(std::stod(Fields[0])*D2R)));
        }
    }
    WAfile.close();

    // Specify compiler flags. These are read in from a file named "SolverParams.bemin"
    std::vector<Real> EnvVars;
    FilePath = std::string("Input") + '/' + std::string("SolverParams.bemin");
    std::ifstream SPfile;
    SPfile.open(FilePath);
    if (SPfile.is_open())
    {
        while ( std::getline(SPfile,line) )
        {
            std::vector<std::string> Fields = Split(line,' ');              // Split line into segments
            if (Fields[1] == "Irregular"){
                if (Fields[0] == "TRUE")    Params.push_back(BEMUse::Parameter("IFR",true));
            }
            if (Fields[1] == "Depth"){
                if (Fields[0] == "INF")     {}// Do nothing (Automatically set to inifite depth)}
                else                        std::cout << "BEMUse is not yet configured for finite depth problems. Analysis shall continue assuming infinite depth.\n";
            }
            if (Fields[1] == "Kochin")          Params.push_back(BEMUse::Parameter("NKochin",std::stoi(Fields[0])));
            if (Fields[1] == "Density")         Params.push_back(BEMUse::Parameter("Density",std::stoi(Fields[0])));
            if (Fields[1] == "Acceleration")    Params.push_back(BEMUse::Parameter("Gravity",Real(std::stod(Fields[0]))));

        }
    }
    SPfile.close();

    //--- Create solver object
    Solver = new BEMUse::Hydrodynamic_Radiation_Solver();

    //--- Set solver parameters
    Solver->Set_Parameters(Params);
    std::string OutputPath = argv[2];
    Solver->Set_OutputFilePath(OutputPath);

    //--- Initialize solver
    Solver->Setup(Boundary);

    std::cout << "Setup of Hydrodynamic radiation solver complete.\n";
}

//--- Execute

void BEMUse_Console::Execute()
{
    // If the solver has not been created, no need to proceed
    if (Solver==NULL){
        if (Boundary!=NULL) delete Boundary;
        std::cout << "Solver undefined. Analysis aborted. Geometry objects cleared.\n";
        return;
    }

    // The sim has been defined, now simply execute
    for (int i=0; i<Frequencies.size(); i++)
    {
        Solver->Set_Real(Frequencies[i]);
        Solver->Solve();
        std::cout << "Frequency " << i+1 << " of " << Frequencies.size() << " complete.\n";
    }

    delete Solver;
    delete Boundary;
    std::cout << "Analysis complete. Geometry objects cleared.\n";
}
