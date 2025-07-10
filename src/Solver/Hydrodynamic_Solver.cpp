//---------------------------------------------------------------
//------------------Hydro solver Functions-----------------------
//---------------------------------------------------------------

#include "Hydrodynamic_Solver.h"

namespace BEMUse
{

//--- Problem setup
void Hydrodynamic_Radiation_Solver::Create_Panels(Boundary *B)
{
    // The panels for the analysis are generated here. For hydrodynamic problems we need to use mirrored panels.

    //--- Generate mirrored panels
    B->Generate_Symmetry_Nodes(ZPlane);
    B->Generate_Symmetry_Elements();

    //--- Extract panels
    std::vector<SP_Geo> Geo, GeoRefl;
    B->Get_Elements(Geo);
    B->Get_Symm_Elements(GeoRefl);
    B->Get_FreeSurface_Elements(FS_Geo);
    NPA = Geo.size();               // Set number of geometry surface panels here

    //--- Are we carrying out irregular freqency removal? Then include also the auxiliary panels & reflected auxiliary panels
    if (IFR)
    {
        //--- Generate mirrored auxiliary panels
        B->Generate_Symmetry_Aux_Nodes(ZPlane);
        B->Generate_Symmetry_Aux_Elements();

        //--- Extract auxiliary mirrored panels
        B->Get_Aux_Elements(Geo);
        B->Get_Aux_Symm_Elements(GeoRefl);
    }

    for (int i=0; i<Geo.size(); i++){
        if (Geo[i]->Get_N()==3)     Source_Panels.push_back(std::make_shared<FlatSourceTriPanel>(Geo[i]));
        if (Geo[i]->Get_N()==4)     Source_Panels.push_back(std::make_shared<FlatSourceQuadPanel>(Geo[i]));

        if (Geo[i]->Get_N()==3)     Refl_Source_Panels.push_back(std::make_shared<FlatSourceTriPanel>(GeoRefl[i]));
        if (Geo[i]->Get_N()==4)     Refl_Source_Panels.push_back(std::make_shared<FlatSourceQuadPanel>(GeoRefl[i]));

        if (Geo[i]->Get_N()==3)     Wave_Panels.push_back(std::make_shared<FlatWaveTriPanel>(GeoRefl[i]));
        if (Geo[i]->Get_N()==4)     Wave_Panels.push_back(std::make_shared<FlatWaveQuadPanel>(GeoRefl[i]));
    }

    NPTot = Source_Panels.size();

    // For simplicity I will simply store the source panels, as these are where the BC is calculated
    StdAppend(Panels,Source_Panels);

    // Now specify panel nodes (for BC later)
    for (SP_Panel P : Panels)   Panel_Nodes.push_back(P->Get_Geo()->Centroid);

//    std::cout   << "NSource panels = " <<  Source_Panels.size()
//                << " NWave Symmetry panels = " << Refl_Source_Panels.size()
//                << " N Wave panels = " << Wave_Panels.size()
//                << " NPA = " << NPA << " NPTot = "  << NPTot << std::endl;

}

void Hydrodynamic_Radiation_Solver::Prepare_Linear_System()
{
    // Prepare the linear system

    SSrcMat =       CMatrix::Zero(NPTot,NPTot);
    SSrcReflMat =   CMatrix::Zero(NPTot,NPTot);
    DSrcMat =       CMatrix::Zero(NPTot,NPTot);
    DSrcReflMat =   CMatrix::Zero(NPTot,NPTot);
    ReflNormMat =   CMatrix::Zero(NPTot,1);

    OpenMPfor
    for (int S=0; S<NPTot; S++){

        Vector3 PanNorm = Source_Panels[S]->Centroid()->Z_Axis_Global();
        ReflNormMat(S) = PanNorm(2)*2.0;

        for (int R=0; R<NPTot; R++)
        {
            // Source terms
            Real s,d,sr,dr;
            Source_Panels[S]->Inf_SingleDoubleLayer(BC_Pos[R],s,d);         // Rankine source
            Refl_Source_Panels[S]->Inf_SingleDoubleLayer(BC_Pos[R],sr,dr);  // Reflected rankine source

            SSrcMat(R,S) = CReal(-s,0);
            SSrcReflMat(R,S) = CReal(-sr,0);
            DSrcMat(R,S) = CReal(d,0);
            DSrcReflMat(R,S) = CReal(dr,0);
        }
    }

//    for (int i=0; i<NPTot; i++) std::cout << SSrcMat(i,i).real() << " " << SSrcMat(i,i).imag() <<" "<< DSrcMat(i,i).real() << " " << DSrcMat(i,i).imag() << std::endl;

}

void Hydrodynamic_Radiation_Solver::Set_RHS_Mat()
{
    // This sets the BC for the problem.

    N_k_Mat = CMatrix::Zero(NPTot,NDOF);    // Reset

    for (int i=0; i<NPTot; i++)
    {
        Vector3 P = Panel_Nodes[i]->Position_Global();
        Vector3 N = Panel_Nodes[i]->Z_Axis_Global();
        Vector3 PCN = P.cross(N);

        N_k_Mat(i,0) = N(0);
        N_k_Mat(i,1) = N(1);
        N_k_Mat(i,2) = N(2);
        N_k_Mat(i,3) = PCN(0);
        N_k_Mat(i,4) = PCN(1);
        N_k_Mat(i,5) = PCN(2);
    }

    DPhi_J_DN = N_k_Mat.block(0,0,NPA,NDOF);    // Gradient of radiation potential
}

void Hydrodynamic_Radiation_Solver::Prepare_PostProcessing_Mats()
{
    // The surface integrals used to calcualte the solutions are able to be expressed much more simply as large matrix multiplications
    // This is done here to improve overview of the code

    AddedMassMat = CMatrix::Zero(NDOF,NPA);
    PanArea = Matrix::Zero(NPA,1);

    OpenMPfor
    for (int i=0; i<NPA; i++)
    {
        Vector3 P = Panel_Nodes[i]->Position_Global();
        Vector3 N = Panel_Nodes[i]->Z_Axis_Global();
        Vector3 PCN = P.cross(N);
        Real dS = Panels[i]->Get_Geo()->Area;

        PanArea(i) = dS;
        AddedMassMat(0,i) = N(0)*dS;
        AddedMassMat(1,i) = N(1)*dS;
        AddedMassMat(2,i) = N(2)*dS;
        AddedMassMat(3,i) = PCN(0)*dS;
        AddedMassMat(4,i) = PCN(1)*dS;
        AddedMassMat(5,i) = PCN(2)*dS;
    }

}

void Hydrodynamic_Radiation_Solver::Prepare_Linear_System_Wave_Terms()
{
    // Prepare the linear system for the wave terms.
    // These coefficients must be calculated at each frequency.

    // I will jump out if we are inspecting inf or zero frequencies
    if (Omega==0.0)     return;
    if (Omega==FreqInf) return;

    // Set matrices to zero
    SWaveMat = CMatrix::Zero(NPTot,NPTot);
    DWaveMat = CMatrix::Zero(NPTot,NPTot);

    //Set Wavenumber parameter in panels
    for (SP_Panel P : Wave_Panels)      P->k = Kappa;

    OpenMPfor
    for (int S=0; S<NPTot; S++){
        for (int R=0; R<NPTot; R++)
        {
            CReal ws,wd;
            Wave_Panels[S]->CInf_SingleDoubleLayerQuad(BC_Pos[R],ws,wd);    // Wave term
            SWaveMat(R,S) = -ws;
            DWaveMat(R,S) = wd + SSrcReflMat(R,S)*ReflNormMat(S)*Kappa;
        }
    }
}

void Hydrodynamic_Radiation_Solver::Prepare_FS_Linear_System(Boundary *B)
{
    // The calculation of the external wave system is accomplished by determining the value of the potential
    // at the external free surface positions. This can be precalculated as with the surface problem solution.

    // Generate list of node positions
    B->Get_FreeSurface_Nodes(Wave_Nodes);
    if (Wave_Nodes.empty()) return;

    int NWN = Wave_Nodes.size();
    std::vector<Vector3> WaveNodePos;
    for (SP_Node N: Wave_Nodes) WaveNodePos.push_back(N->Position_Global());

    // Create Array
    SMatExt = CMatrix::Zero(NWN,NPA);
    SMatReflExt = CMatrix::Zero(NWN,NPA);
    DMatExt = CMatrix::Zero(NWN,NPA);
    DMatReflExt = CMatrix::Zero(NWN,NPA);

    OpenMPfor
        for (int S=0; S<NPA; S++){
        for (int R=0; R<NWN; R++)
        {
            // Source terms
            Real s,d,sr,dr;
            Source_Panels[S]->Inf_SingleDoubleLayer(WaveNodePos[R],s,d);         // Rankine source
            Refl_Source_Panels[S]->Inf_SingleDoubleLayer(WaveNodePos[R],sr,dr);  // Reflected rankine source
            SMatExt(R,S) = CReal(-s,0);
            SMatReflExt(R,S) = CReal(-sr,0);
            DMatExt(R,S) = CReal(d,0);
            DMatReflExt(R,S) = CReal(dr,0);
        }
    }
}

void Hydrodynamic_Radiation_Solver::Prepare_FS_Linear_System_Wave_Terms()
{
    // Prepare the frequency-dependent coefficients of the linear system for the
    // calculation of the potential on the free surface.

    // Generate list of node positions
    if (Wave_Nodes.empty()) return;             // No external nodes, ignore!

    int NWN = Wave_Nodes.size();
    std::vector<Vector3> WaveNodePos;
    for (SP_Node N: Wave_Nodes) WaveNodePos.push_back(N->Position_Global());

    // Create wave term array
    SWaveMatExt = CMatrix::Zero(NWN,NPA);
    DWaveMatExt = CMatrix::Zero(NWN,NPA);

    OpenMPfor
    for (int S=0; S<NPA; S++){
        for (int R=0; R<NWN; R++)
        {
            CReal ws,wd;
            Wave_Panels[S]->CInf_SingleDoubleLayerQuad(WaveNodePos[R],ws,wd);    // Wave term
            SWaveMatExt(R,S) = -ws;
            DWaveMatExt(R,S) = wd + SMatReflExt(R,S)*ReflNormMat(S)*Kappa;
        }
    }
}

//--- Setup
void Hydrodynamic_Radiation_Solver::Setup(Boundary *B)
{
    // Prepare system
    Create_Panels(B);
    Specify_BC_Const_Pans(B);
    Prepare_Linear_System();

    // Prepare matrices for post processing
    Set_RHS_Mat();                      // Set BCs for each DOF
    Prepare_PostProcessing_Mats();

    // Prepare exterior free surface
    Prepare_FS_Linear_System(B);

    // Specify storage arrays
    NBeta = BetaArray.size();
    FK_i = CMatrix::Zero(NDOF,NBeta);
    F_i = CMatrix::Zero(NDOF,NBeta);
    RAO_i = CMatrix::Zero(NDOF,NBeta);

    // Generate initial output file
    Calculate_Mass_Matrix(B);
    Calculate_Hydrostatic_Stiffness_Matrix(B);
    Generate_Output_File(B);

    // Solve problem for very large frequency (necessary for added mass in Cummins equation)
    Frequency = FreqInf;
    Solve();
    Update_Output_File();
}

void Hydrodynamic_Radiation_Solver::Calculate_Mass_Matrix(Boundary *B)
{
    // This calculates the mass matrix based on the values specified for the geometry

    // First carry out a check here to see if the mass matrix was previously set
    if (!B->Mass_Mat.isZero())
    {
        M_ij = B->Mass_Mat;
        return;
    }

    if (B->Mass==0) B->Mass = B->Volume*Rho_wat;            // Assume water density

    // Upper left
    M_ij(0,0) = B->Mass;
    M_ij(1,1) = B->Mass;
    M_ij(2,2) = B->Mass;

    // Upper right
    M_ij(0,4) = B->Mass*B->Centre_Gravity(2);
    M_ij(1,3) = -B->Mass*B->Centre_Gravity(2);
    M_ij(0,5) = -B->Mass*B->Centre_Gravity(1);
    M_ij(2,3) = B->Mass*B->Centre_Gravity(1);
    M_ij(1,5) = B->Mass*B->Centre_Gravity(0);
    M_ij(2,4) = -B->Mass*B->Centre_Gravity(0);

    // Lower left
    M_ij.block(3,0,3,3) = -M_ij.block(0,3,3,3);

    // Lower right
    M_ij(3,3) = B->Mass*B->Radii_Gyration(0,0)*fabs(B->Radii_Gyration(0,0));
    M_ij(3,4) = B->Mass*B->Radii_Gyration(0,1)*fabs(B->Radii_Gyration(0,1));
    M_ij(3,5) = B->Mass*B->Radii_Gyration(0,2)*fabs(B->Radii_Gyration(0,2));
    M_ij(4,3) = B->Mass*B->Radii_Gyration(1,0)*fabs(B->Radii_Gyration(1,0));
    M_ij(4,4) = B->Mass*B->Radii_Gyration(1,1)*fabs(B->Radii_Gyration(1,1));
    M_ij(4,5) = B->Mass*B->Radii_Gyration(1,2)*fabs(B->Radii_Gyration(1,2));
    M_ij(5,3) = B->Mass*B->Radii_Gyration(2,0)*fabs(B->Radii_Gyration(2,0));
    M_ij(5,4) = B->Mass*B->Radii_Gyration(2,1)*fabs(B->Radii_Gyration(2,1));
    M_ij(5,5) = B->Mass*B->Radii_Gyration(2,2)*fabs(B->Radii_Gyration(2,2));
}

void Hydrodynamic_Radiation_Solver::Calculate_Hydrostatic_Stiffness_Matrix(Boundary *B)
{
    // We shall assume Centre of flotation is at zero

    // First carry out a check here to see if the stiffness matrix was previously set
    if (!B->Stiffness_Mat.isZero())
    {
        C_ij = B->Stiffness_Mat;
        return;
    }

    //--- Calculate the moments of the free surface area
    std::vector<SP_Geo> Aux_Els;
    B->Get_Aux_Elements(Aux_Els);
    Real AFS = 0, MOAx = 0, MOAy = 0, MOAxx = 0, MOAxy = 0, MOAyx = 0, MOAyy = 0;
    for (SP_Geo G : Aux_Els){
        Vector3 Pos = G->Centroid->Position_Global();
        Real A = G->Area;

        AFS += A;
        MOAy += A*Pos(1);
        MOAx += A*Pos(0);
        MOAxy += A*Pos(0)*Pos(1);
        MOAyx += A*Pos(1)*Pos(0);
        MOAxx += A*Pos(0)*Pos(0);
        MOAyy += A*Pos(1)*Pos(1);
    }

    //--- This is taken from the WAMIT user manual

    Real    g = Grav;
    Real    rho = Rho_wat;
    Real    M = B->Mass;
    Real    V = B->Volume;

    Real xb = B->Centre_Buoyancy(0);
    Real yb = B->Centre_Buoyancy(1);
    Real zb = B->Centre_Buoyancy(2);

    Real xg = B->Centre_Gravity(0);
    Real yg = B->Centre_Gravity(1);
    Real zg = B->Centre_Gravity(2);

    C_ij(2,2) =  Rho*g*AFS;
    C_ij(2,3) =  Rho*g*MOAy;
    C_ij(2,4) = -Rho*g*MOAx;
    C_ij(3,3) =  Rho*g*MOAyy + rho*g*V*zb - M*g*zg;
    C_ij(3,4) = -Rho*g*MOAxy;
    C_ij(3,5) = -Rho*g*V*xb + M*g*xg;
    C_ij(4,4) =  Rho*g*MOAxx + rho*g*V*zb - M*g*zg;
    C_ij(4,5) = -Rho*g*V*yb + M*g*yg;

    // Now apply symmetries
    C_ij(2,3) = C_ij(3,2);
    C_ij(2,4) = C_ij(4,2);
    C_ij(3,4) = C_ij(4,3);
}

//--- Solution
void Hydrodynamic_Radiation_Solver::Solve()
{
    // This solve function calculates the hydrodynamic coefficients for a single frequency

    //--- Set simulation parameters
    Omega = Frequency*TwoPIinv;
    if (Omega==0 && IFR)    Omega = 1.0e-3; // Avoid erroneous values for small arguments

    Period = 1.0/Frequency;
    if (Depth == Infinite)   Kappa = Omega*Omega/Gravity;
    if (Depth == Finite)     Kappa = Calc_Root_Finite_Depth(Omega, H);

    //--- Compose full influence coefficient matrices
    Prepare_Linear_System_Wave_Terms();
    CMatrix SMat, DMat;
    if (Omega==0.0)                         // Case 1: Frequency zero
    {
        SMat = SSrcMat + SSrcReflMat;
        DMat = DSrcMat + DSrcReflMat;
    }
    else if (Omega==FreqInf)                // Case 2: Frequency infinite (for added mass terms)
    {
        SMat = SSrcMat - SSrcReflMat;
        DMat = DSrcMat - DSrcReflMat;
    }
    else                                    // Case 3: Full wave effects
    {
        SMat = SSrcMat + SSrcReflMat + SWaveMat;
        DMat = DSrcMat + DSrcReflMat + DWaveMat;
    }
    PPLU.compute(DMat);                     // Prepare linear solver for solution
//    CG.compute(DMat);                     // Prepare linear solver for solution
//    DGMRES.compute(DMat);
//    GMRES.compute(DMat);
//    VisMat = DMat;

    //--- Calc radiation solution   (one solution per DOF)
    CMatrix RHSTemp = SMat*N_k_Mat;      // Source matrix term
    CMatrix Phi_J = PPLU.solve(RHSTemp);                       // Solution using a partial piv Lu decomposition
//    CMatrix Phi_J = CG.solve(RHSTemp);                       // Solution using a partial piv Lu decomposition
//    CMatrix Phi_J = DGMRES.solve(RHSTemp);
//    CMatrix Phi_J = GMRES.solve(RHSTemp);

    //--- Calc diffraction solution (one solution per incoming wave angle Beta)
    Set_Incident_Potential_Mats();
    CMatrix RHSDiff = SMat*(-DPhi_I_DN.cwiseProduct(PanArea));
    CMatrix Phi_S = PPLU.solve(RHSDiff);
//    CMatrix Phi_S = CG.solve(RHSDiff);
//    CMatrix Phi_S = DGMRES.solve(RHSDiff);
//    CMatrix Phi_S = GMRES.solve(RHSDiff);
    CMatrix Phi_D = Phi_S+Phi_I;

    //--- Discard solutions which are not of interest (only panels on the surface)
    Phi_J.conservativeResize(NPA,NDOF);         // Radiation potential
    Phi_S.conservativeResize(NPA,NBeta);        // Scatter potential
    Phi_D.conservativeResize(NPA,NBeta);        // Diffraction potential
    Phi_I.conservativeResize(NPA,NBeta);        // Incident potential
    DPhi_I_DN.conservativeResize(NPA,NBeta);    // Gradient of incident potential

    //--- Post processing using large matrix operations

    // Added mass
    CMatrix AM_Mat = AddedMassMat*Phi_J;
    A_ij = AM_Mat.real();
    B_ij = AM_Mat.imag();

    //--- Set a flag... If the infinite frequency analysis is being carried out, jump out here.
    //--- We only need radiation force at high frequency for the Cummins equation.
    if (Frequency==FreqInf) return;

    // Excitation forces
    CMatrix F1 = AddedMassMat*Phi_I;
    CMatrix F2 = Phi_J.transpose()*(DPhi_I_DN.cwiseProduct(PanArea));
    FK_i = -Im*Omega*F1;                // Froude-Krylov forces
    SC_i =  Im*Omega*F2;                // Scattering forces
    CMatrix FSup = FK_i + SC_i;
    F_i = FSup.conjugate();             // Excitation forces

    //-- Correct in case of zero frequency
    if (Omega==0.0){
        F_i = CMatrix::Zero(6,NBeta);
        FK_i = CMatrix::Zero(6,NBeta);
    }

    // Response amplitude operator
    CMatrix LHS = -Omega*Omega*(A_ij+M_ij/Rho_wat) + Omega*Im*B_ijOm + C_ij/Rho_wat;
    RAO_i = (LHS.inverse())*F_i;

    //--- Kochin integral (all notation follows Newman 1967)

    //------------------- Newman approach
    Set_Kochin_Mats();
    KochinRad = KochinE*DPhi_J_DN - KochindEdn*Phi_J;           // Kochin radiation [NKoch,NDOF]
    KochinDiff = KochinE*DPhi_I_DN - KochindEdn*Phi_S;          // Kochin diffraction [NKoch,NBeta]
    // These terms must be integrated azimuthally in order to calculate mean drift terms.
    // This is left to the user as the RAOs must be used.

    //--- Calc free surface elevation
    if (!Wave_Nodes.empty()){
        Prepare_FS_Linear_System_Wave_Terms();
        CMatrix WSMat = SMatExt + SMatReflExt + SWaveMatExt;
        CMatrix WDMat = DMatExt + DMatReflExt + DWaveMatExt;
        ExtRadMat = -Im*Omega/Gravity*(WSMat*DPhi_J_DN - WDMat*Phi_J);
        ExtDiffMat = -Im*Omega/Gravity*(WSMat*DPhi_I_DN + WDMat*Phi_S);
    }

    //--- Update output files
    Update_Output_File();
    if (!Wave_Nodes.empty()) Export_Wave_Height();

    //--- Store results for visualisation
    RadSolArray.push_back(Phi_J);
    DiffSolArray.push_back(Phi_S);
    FS_Rad_SolArray.push_back(ExtRadMat);
    FS_Scat_SolArray.push_back(ExtDiffMat);
}

Real Hydrodynamic_Radiation_Solver::Calc_Root_Finite_Depth(Real &O, Real &H)
{
    // This function makes use of formulas (25) and (27) from:
    // Newman, J. N. 1990 “Numerical solutions of the water-wave dispersion relation,” Applied Ocean Research, 12, 14-18.
    // In order to calculate the roots of the equation κ tanh κh = ω^2/g

    Real x = O*O*H/Gravity;
//    Real An[9] = {0.03355458, 0.03262249, -0.00088239, 0.00004620, -0.00000303, 0.00000034, -0.00000007, 0.00000003, -0.00000001};
    Real Bn[6] = {0.000000122, 0.073250017, -0.009899981, 0.002640863, 0.000829239, -0.000176411};
    Real Cn[9] = {1.0, -0.33333372, -0.01109668, 0.01726435, 0.01325580, -0.00116594, 0.00829006, -0.01252603, 0.00404923};

    Real y = 0.0;
    if (x <= 2.0){
        Real Den = 0.0;
        for (int i=0; i<=8; i++) Den += Cn[i]*pow( 0.5*x , i);
        y = sqrt(x)/Den;                                                        // Maximum error e = 1.0 x 10^-8.
    }
    else {
        y += x;
        for (int i=0; i<=5; i++) y += Bn[i]*pow( 0.5*x*exp(4.0-2.0*x) , i);     // Maximum error e = 1.2 x 10^-7.
    }
    return y/H;

    std::cout << "Hydrodynamic_Boundary::Calc_Root_Finite_Depth: Solution of Kappa for finite depth problems must be implemented.\n";
}

void Hydrodynamic_Radiation_Solver::Calc_Phi_Incident(Real &B, Vector3 &P_Glob, CReal &Phi_i, CReal &dPhi_iX,  CReal &dPhi_iY,  CReal &dPhi_iZ)
{
    // This function calcualtes the value of an incident wave based on the wave
    // parameters and the type of model
    Real x = P_Glob(0), y = P_Glob(1), z = P_Glob(2);
    Real PS = 0;
    Real CB = cos(B+PS), SB = sin(B+PS);
    CReal Phi_I = Im*Gravity/Omega*exp(Im*Kappa*(x*CB+y*SB));

    if (Depth == Infinite)
    {
        // See equation 2.1 WAMIT theory manual
        Real Z = exp(Kappa*z);
        Phi_i = Z*Phi_I;
        dPhi_iX = Im*Kappa*CB*Phi_i;
        dPhi_iY = Im*Kappa*SB*Phi_i;
        dPhi_iZ = Kappa*Phi_i;
    }
    if (Depth == Finite)
    {
        // See equation 2.3 WAMIT theory manual
        Real Z = cosh(Kappa*(z+H))/cosh(Kappa*H);
        Real dZ = Kappa*sinh(Kappa*(z+H))/cosh(Kappa*H);

        Phi_i = Im*Z*Phi_I;
        dPhi_iX = Im*Kappa*CB*Phi_i;
        dPhi_iY = Im*Kappa*SB*Phi_i;
        dPhi_iZ = dZ*Phi_i;
    }
}

void Hydrodynamic_Radiation_Solver::Set_Incident_Potential_Mats()
{
    // The incident potential & its gradient on the surface must be calculated.
    Phi_I = CMatrix::Zero(NPTot,NBeta);
    DPhi_I_DN = CMatrix::Zero(NPTot,NBeta);

    OpenMPfor
    for (int P=0; P<NPTot; P++){              // Loop panels

        Vector3 Pos = Panel_Nodes[P]->Position_Global();     // Node position
        Vector3 Norm = Panel_Nodes[P]->Z_Axis_Global();            // Norm of the panel

        for (int B = 0; B<NBeta; B++){      // Loop Incoming Wave Angles

            CReal pi, dpx, dpy, dpz;
            Calc_Phi_Incident(BetaArray[B], Pos, pi, dpx, dpy, dpz);

            Phi_I(P,B) = pi;
            DPhi_I_DN(P,B) = Norm(0)*dpx + Norm(1)*dpy + Norm(2)*dpz;
        }
    }
}

void Hydrodynamic_Radiation_Solver::Set_Kochin_Mats()
{
    // These store the E and dE/dN values from the Kochin function which are used to calculate the Kochin functions

    // Allocate arrays
    KochinE = CMatrix::Zero(NKoch,NPA);
    KochindEdn = CMatrix::Zero(NKoch,NPA);

    KochinPiE = CMatrix::Zero(1,NPA);
    KochinPidEdn = CMatrix::Zero(1,NPA);

    // Create list of surface positions, normals
    std::vector<Vector3> Pos, Norm;
    for (int i=0; i<NPA; i++){
        Pos.push_back(Panel_Nodes[i]->Position_Global());
        Norm.push_back(Panel_Nodes[i]->Z_Axis_Global());
    }

    OpenMPfor
    for (int i=0; i<NKoch+1; i++)
    {
        Real Theta = (0.5+1.0*i)*TwoPI/NKoch;
        if (i==NKoch) Theta = Beta + PI;                   // Constant theta term (see Newman)
        Real CTH = cos(Theta), STH = sin(Theta);

        for (int P=0; P<NPA; P++)
        {
            Real x = Pos[P](0), y = Pos[P](1), z = Pos[P](2);
            CReal(Kappa*z, Kappa*(x*CTH + y*STH));
            CReal E = exp(E);
            CReal dEdx = Im*Kappa*CTH*E;
            CReal dEdy = Im*Kappa*STH*E;
            CReal dEdz = Kappa*E;

            if (i==NKoch){
                KochinPiE(P) =  E*PanArea(P);
                KochinPidEdn(P) = (Norm[P](0)*dEdx + Norm[P](1)*dEdy + Norm[P](2)*dEdz)*PanArea(P);
            }
            else{
                KochinE(i,P) = E*PanArea(P);
                KochindEdn(i,P) = (Norm[P](0)*dEdx + Norm[P](1)*dEdy + Norm[P](2)*dEdz)*PanArea(P);
            }
        }
    }
}

}
