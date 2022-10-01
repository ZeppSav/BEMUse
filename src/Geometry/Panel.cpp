//-----------------------------------------------------------------------------
//-------------------------Panel Element Functions-----------------------------
//-----------------------------------------------------------------------------

#include "Panel.h"
#include "../Kernels.cpp"

namespace BEMUse
{

//-----------------------
// Flat Tri Panel
//-----------------------

// The source kernel appears so often, I will just code it into the parent class to avoid repetition

void FlatTriPanel::Inf_SourceKernel(const Vector3 &P, Real &S, Real &D)
{
    // Analytical expression for source panel.

    Real m[3], e[3], h[3], r[3];

    // Calculate geometric parameters.

    Vector3 CLoc = Centroid()->Position_Global();
    Matrix3 OLoc = Centroid()->OrientMat_Global();

    Vector3 P0 = OLoc*(Geo->Nodes[0]->Position_Global() - CLoc);
    Vector3 P1 = OLoc*(Geo->Nodes[1]->Position_Global() - CLoc);
    Vector3 P2 = OLoc*(Geo->Nodes[2]->Position_Global() - CLoc);
    Vector3 PL = OLoc*(P-CLoc);

    m[0] = (P1(1)-P0(1))/(P1(0)-P0(0));
    m[1] = (P2(1)-P1(1))/(P2(0)-P1(0));
    m[2] = (P0(1)-P2(1))/(P0(0)-P2(0));

    Real DX0 = PL(0)-P0(0), DY0 = PL(1)-P0(1);
    Real DX1 = PL(0)-P1(0), DY1 = PL(1)-P1(1);
    Real DX2 = PL(0)-P2(0), DY2=  PL(1)-P2(1);

    Real Z2 = PL(2)*PL(2);

    e[0] =  DX0*DX0 + Z2;
    e[1] =  DX1*DX1 + Z2;
    e[2] =  DX2*DX2 + Z2;

    h[0] = DX0*DY0;
    h[1] = DX1*DY1;
    h[2] = DX2*DY2;

    r[0] = sqrt( DX0*DX0 + DY0*DY0 + Z2 );
    r[1] = sqrt( DX1*DX1 + DY1*DY1 + Z2 );
    r[2] = sqrt( DX2*DX2 + DY2*DY2 + Z2 );

    Real PhiDip = 0;

    PhiDip += atan( (m[0]*e[0]-h[0]) / (PL(2)*r[0]) ) - atan( (m[0]*e[1]-h[1]) / (PL(2)*r[1]) );
    PhiDip += atan( (m[1]*e[1]-h[1]) / (PL(2)*r[1]) ) - atan( (m[1]*e[2]-h[2]) / (PL(2)*r[2]) );
    PhiDip += atan( (m[2]*e[2]-h[2]) / (PL(2)*r[2]) ) - atan( (m[2]*e[0]-h[0]) / (PL(2)*r[0]) );

    Real d[3], B[3], L[3];

    d[0] = sqrt( (P1(0)-P0(0))*(P1(0)-P0(0)) + (P1(1)-P0(1))*(P1(1)-P0(1)) );
    d[1] = sqrt( (P2(0)-P1(0))*(P2(0)-P1(0)) + (P2(1)-P1(1))*(P2(1)-P1(1)) );
    d[2] = sqrt( (P0(0)-P2(0))*(P0(0)-P2(0)) + (P0(1)-P2(1))*(P0(1)-P2(1)) );

    B[0] = ( DX0*(P1(1)-P0(1))-DY0*(P1(0)-P0(0)) )/d[0];
    B[1] = ( DX1*(P2(1)-P1(1))-DY1*(P2(0)-P1(0)) )/d[1];
    B[2] = ( DX2*(P0(1)-P2(1))-DY2*(P0(0)-P2(0)) )/d[2];

    L[0] = log((r[0]+r[1]+d[0]) / (r[0]+r[1]-d[0]));
    L[1] = log((r[1]+r[2]+d[1]) / (r[1]+r[2]-d[1]));
    L[2] = log((r[2]+r[0]+d[2]) / (r[2]+r[0]-d[2]));

    // Set outputs
    D = PhiDip*FourPIinv;
    S = (PL(2)*PhiDip - (B[0]*L[0] + B[1]*L[1] + B[2]*L[2]))*FourPIinv;
}

//-----------------------
// Flat Quad Panel
//-----------------------

// The source kernel appears so often, I will just code it into the parent class to avoid repetition

void FlatQuadPanel::Inf_SourceKernel(const Vector3 &P, Real &S, Real &D)
{
    // Analytical function.
    Real m[4], e[4], h[4], r[4];

    // Calculate geometric parameters.

    Vector3 CLoc = Centroid()->Position_Global();
    Matrix3 OLoc = Centroid()->OrientMat_Global();

    Vector3 P0 = OLoc*(Geo->Nodes[0]->Position_Global() - CLoc);
    Vector3 P1 = OLoc*(Geo->Nodes[1]->Position_Global() - CLoc);
    Vector3 P2 = OLoc*(Geo->Nodes[2]->Position_Global() - CLoc);
    Vector3 P3 = OLoc*(Geo->Nodes[3]->Position_Global() - CLoc);
    Vector3 PL = OLoc*(P-CLoc);

    m[0] = (P1(1)-P0(1))/(P1(0)-P0(0));
    m[1] = (P2(1)-P1(1))/(P2(0)-P1(0));
    m[2] = (P3(1)-P2(1))/(P3(0)-P2(0));
    m[3] = (P0(1)-P3(1))/(P0(0)-P3(0));

    Real DX0 = PL(0)-P0(0), DY0 = PL(1)-P0(1);
    Real DX1 = PL(0)-P1(0), DY1 = PL(1)-P1(1);
    Real DX2 = PL(0)-P2(0), DY2=  PL(1)-P2(1);
    Real DX3 = PL(0)-P3(0), DY3 = PL(1)-P3(1);

    Real Z2 = PL(2)*PL(2);

    e[0] =  DX0*DX0 + Z2;
    e[1] =  DX1*DX1 + Z2;
    e[2] =  DX2*DX2 + Z2;
    e[3] =  DX3*DX3 + Z2;

    h[0] = DX0*DY0;
    h[1] = DX1*DY1;
    h[2] = DX2*DY2;
    h[3] = DX3*DY3;

    r[0] = sqrt( DX0*DX0 + DY0*DY0 + Z2 );
    r[1] = sqrt( DX1*DX1 + DY1*DY1 + Z2 );
    r[2] = sqrt( DX2*DX2 + DY2*DY2 + Z2 );
    r[3] = sqrt( DX3*DX3 + DY3*DY3 + Z2 );

    Real PhiDip = 0;

    PhiDip += atan( (m[0]*e[0]-h[0]) / (PL(2)*r[0]) ) - atan( (m[0]*e[1]-h[1]) / (PL(2)*r[1]) );
    PhiDip += atan( (m[1]*e[1]-h[1]) / (PL(2)*r[1]) ) - atan( (m[1]*e[2]-h[2]) / (PL(2)*r[2]) );
    PhiDip += atan( (m[2]*e[2]-h[2]) / (PL(2)*r[2]) ) - atan( (m[2]*e[3]-h[3]) / (PL(2)*r[3]) );
    PhiDip += atan( (m[3]*e[3]-h[3]) / (PL(2)*r[3]) ) - atan( (m[3]*e[0]-h[0]) / (PL(2)*r[0]) );

    Real d[4], B[4], L[4];

    d[0] = sqrt( (P1(0)-P0(0))*(P1(0)-P0(0)) + (P1(1)-P0(1))*(P1(1)-P0(1)) );
    d[1] = sqrt( (P2(0)-P1(0))*(P2(0)-P1(0)) + (P2(1)-P1(1))*(P2(1)-P1(1)) );
    d[2] = sqrt( (P3(0)-P2(0))*(P3(0)-P2(0)) + (P3(1)-P2(1))*(P3(1)-P2(1)) );
    d[3] = sqrt( (P0(0)-P3(0))*(P0(0)-P3(0)) + (P0(1)-P3(1))*(P0(1)-P3(1)) );

    B[0] = ( DX0*(P1(1)-P0(1))-DY0*(P1(0)-P0(0)) )/d[0];
    B[1] = ( DX1*(P2(1)-P1(1))-DY1*(P2(0)-P1(0)) )/d[1];
    B[2] = ( DX2*(P3(1)-P2(1))-DY2*(P3(0)-P2(0)) )/d[2];
    B[3] = ( DX3*(P0(1)-P3(1))-DY3*(P0(0)-P3(0)) )/d[3];

    L[0] = log((r[0]+r[1]+d[0]) / (r[0]+r[1]-d[0]));
    L[1] = log((r[1]+r[2]+d[1]) / (r[1]+r[2]-d[1]));
    L[2] = log((r[2]+r[3]+d[2]) / (r[2]+r[3]-d[2]));
    L[3] = log((r[3]+r[0]+d[3]) / (r[3]+r[0]-d[3]));

    // Set outputs
    D = PhiDip*FourPIinv;
    S = (PL(2)*PhiDip - (B[0]*L[0] + B[1]*L[1] + B[2]*L[2] + B[3]*L[3]))*FourPIinv;
}

//-----------------------
// Flat Source Tri Panel
//-----------------------

void FlatSourceTriPanel::Inf_SingleDoubleLayerQuad(const Vector3 &P, Real &S, Real &D)
{
    Real Sout = 0, Dout = 0;
    for (int i=0; i<Geo->QNodes.size(); i++){
        Vector S = Geo->QNodes[i]->Position_Global();
        Vector G = Geo->QNodes[i]->Z_Axis_Global();
        Sout += Geo->QNodes[i]->Weight*Source_Kernel(S,P);
        Dout += Geo->QNodes[i]->Weight*Dipole_Kernel(S,P,G);
    }
    S = Sout;
    D = Dout;
}

//-----------------------
// Flat Source Quad Panel
//-----------------------

void FlatSourceQuadPanel::Inf_SingleDoubleLayerQuad(const Vector3 &P, Real &S, Real &D)
{
    Real Sout = 0, Dout = 0;
    for (int i=0; i<Geo->QNodes.size(); i++){
        Vector S = Geo->QNodes[i]->Position_Global();
        Vector G = Geo->QNodes[i]->Z_Axis_Global();
        Sout += Geo->QNodes[i]->Weight*Source_Kernel(S,P);
        Dout += Geo->QNodes[i]->Weight*Dipole_Kernel(S,P,G);
    }
    S = Sout;
    D = Dout;
}

//-----------------------
// Hydrodynamic panels
//-----------------------

//------Liang/Noblesse Implementation

inline void FreeSurface_Wave_Kernel(const Vector3 &P_glob, const Vector3 &P_Src, const Real &K, CReal &G, CReal &dGX, CReal &dGY, CReal &dGZ)
{
    // This is simply the calculation of the wave component for the infinite depth case between two evaluation points
    // The code has been modified from the original version of Liang & Noblesse for simplicity.
    // It is assumed that the source is at a vertical position z>=0 & the reciever is at a position z<=0

    if (K==0.0)
    {
        G = CReal(0.0,0.0);
        dGX = CReal(0.0,0.0);
        dGY = CReal(0.0,0.0);
        dGZ = CReal(0.0,0.0);
        return;
    }

    Vector3 RVec = P_glob - P_Src;        // Relative position global

    Real XX = RVec(0);
    Real YY = RVec(1);
//    Real V = P_glob(2) + P_Src(2);
    Real V = RVec(2);

    //  Non-dimensionalise input vars
    XX *= K;
    YY *= K;
    V *= K;

//    Kernels::HavelockGF(XX,YY,V,G,dGX,dGY,dGZ);       // Calculate wave component
    HavelockGF(XX,YY,V,G,dGX,dGY,dGZ);       // Calculate wave component

    //  Dimensionalise output vars
    G *= K*FourPIinv;
    dGX *= K*K*FourPIinv;
    dGY *= K*K*FourPIinv;
    dGZ *= K*K*FourPIinv;
}

inline void FreeSurface_L(const Vector3 &P_glob, const Vector3 &P_Src, const Real &K, CReal &G, CReal &dGX, CReal &dGY, CReal &dGZ)
{
    // This is simply the calculation of the wave component for the infinite depth case between two evaluation points

    Vector3 RVec = P_glob - P_Src;        // Relative position global

    Real XX = RVec(0);
    Real YY = RVec(1);
    Real V = P_glob(2) + P_Src(2);

    if (K==0.0)
    {
        G = CReal(0.0,0.0);
        dGX = CReal(0.0,0.0);
        dGY = CReal(0.0,0.0);
        dGZ = CReal(0.0,0.0);
        return;
    }

    //  Non-dimensionalise input vars
    XX *= K;
    YY *= K;
    V *= K;

//    Kernels::HavelockOnlyL(XX,YY,V,G,dGX,dGY,dGZ);       // Calculate wave component

    //  Dimensionalise output vars

    G *= K;
    dGX *= pow(K,2);
    dGY *= pow(K,2);
    dGZ *= pow(K,2);
}

//-----------------------
// Flat Wave Tri Panel
//-----------------------

void  FlatWaveTriPanel::CInf_SingleDoubleLayerQuad(const Vector3 &P, CReal &S, CReal &D)
{
    // Calculate both single and double layer potentials using quadrature rule.

    CReal OS(0,0), OD(0,0);
    for (int i=0; i<Geo->QNodes.size(); i++){
        Vector S = Geo->QNodes[i]->Position_Global();
        Vector N = Geo->QNodes[i]->Z_Axis_Global();
        CReal G, dX, dY, dZ;
        FreeSurface_Wave_Kernel(P,S,k,G,dX,dY,dZ);
//        FreeSurface_L(P,S,k,G,dX,dY,dZ);            // HACK!!! Only for testing
        OS += Geo->QNodes[i]->Weight*G;
        OD += Geo->QNodes[i]->Weight*(dX*N(0)+dY*N(1)+dZ*N(2));
//        OD += Geo->QNodes[i]->Weight*(dX*N(0)+dY*N(1)-dZ*N(2));
    }
    // Set outputs
    S = OS;
    D = OD;
}

//--- Sub panelling (brute force num integration)

void  FlatWaveTriPanel::CInf_SingleDoubleLayerBruteForce(const Vector3 &PGlob, CReal &S, CReal &D)
{
    // Calculate both single and double layer potentials sing brute force subpanneling
    // HACK for panel with normal facing upwards

    int NSP = 100;
    Real DX = 1.0/NSP;

    // Subpanelling (openMP)
    CReal OS[NSP],OD[NSP];
    Real ODA[NSP];
    for (int i=0; i<NSP; i++){
        OS[i] = CReal(0,0);
        OD[i] = CReal(0,0);
        ODA[i] = 0.0;
    }

    Real XL = -0.5;     // Equilateral!!!

    OpenMPfor
    for (int i=0; i<NSP; i++){
        for (int j=0; j<NSP; j++){
            Real xi = (0.5+1.0*i)*DX;
            Real eta = (0.5+1.0*j)*DX;

            Vector3 P(0,0,0);
            P += Geo->Nodes[0]->Position_Global()*(1-xi-eta);
            P += Geo->Nodes[1]->Position_Global()*(eta);
            P += Geo->Nodes[2]->Position_Global()*(xi);
            P(2) = ZP;

            // Scale the area elemnt or skip if in wrong position:
            Real dA = DX*DX;
            if (P(0)<XL-0.001) continue;
            if ((P(0)-XL)<0.5*DX){
                dA *= 0.5;
                if (j==0) dA *= 0.5;    // Top edge     (corner)
                if (i==0) dA *= 0.5;    // Bottom edge  (corner)
            }

            // Surface normals
            Vector3 R_xi = Vector3::Zero() , R_eta = Vector3::Zero();
            R_xi += Geo->Nodes[0]->Position_Global()*(-1);
            R_xi += Geo->Nodes[2]->Position_Global()*(1);
            R_eta += Geo->Nodes[0]->Position_Global()*(-1);
            R_eta += Geo->Nodes[1]->Position_Global()*(1);

            Vector3 R_n = R_xi.cross(R_eta);
            Real Jac = R_n.norm();              // Jacobian of quadrature node. Must be multiplied with weight.

            Vector3 N(0,0,1);
            CReal G, dX, dY, dZ;
            FreeSurface_Wave_Kernel(PGlob,P,k,G,dX,dY,dZ);
//            FreeSurface_L(PGlob,P,k,G,dX,dY,dZ);         // HACK!!! ONLY L Terms!!!
            OS[i] += Jac*G*dA;
            OD[i] += Jac*(dX*N(0)+dY*N(1)+dZ*N(2))*dA;
            ODA[i] += dA;
        }
    }

    CReal OST(0,0),ODT(0,0);
    Real A = 0.0;
    for (int i=0; i<NSP; i++){
        OST += OS[i];
        ODT += OD[i];
        A += ODA[i];
    }
//    qDebug() << "Area check: Should be 0.5 (equilateral), is: " << A;

    // Set outputs
    S = OST;
    D = ODT;
}

//-----------------------
// Flat Wave Quad Panel
//-----------------------

void  FlatWaveQuadPanel::CInf_SingleDoubleLayerQuad(const Vector3 &P, CReal &S, CReal &D)
{
    // Calculate both single and double layer potentials sing quadrature rule.

    CReal OS(0,0), OD(0,0);
    for (int i=0; i<Geo->QNodes.size(); i++){
        Vector S = Geo->QNodes[i]->Position_Global();
        Vector N = Geo->QNodes[i]->Z_Axis_Global();
        CReal G, dX, dY, dZ;
        FreeSurface_Wave_Kernel(P,S,k,G,dX,dY,dZ);
//        FreeSurface_L(P,S,k,G,dX,dY,dZ);         // HACK!!! ONLY L Terms!!!
        OS += Geo->QNodes[i]->Weight*G;
        OD += Geo->QNodes[i]->Weight*(dX*N(0)+dY*N(1)+dZ*N(2));
//        OD += Geo->QNodes[i]->Weight*(dX*N(0)+dY*N(1)-dZ*N(2));
    }
    // Set outputs
    S = OS;
    D = OD;
}

void  FlatWaveQuadPanel::Wave_Kernel_Test()
{
    // THis simply plots as a function of depth and normalized z the wavefunction...
    int NT = 100;
    Real R = 3;
    Vector3 PSrc(0,0,-1);
    for (int i=0; i<NT; i++){

//        Real r = 2.0+3.0*i/(NT-1);
        CReal G1,G2,G3,dX,dY,dZ;
//        FreeSurface_Wave_Kernel(Vector3(r,0,-1),PSrc,k,G1,dX,dY,dZ);
//        FreeSurface_Wave_Kernel(Vector3(r,r,-5),PSrc,k,G2,dX,dY,dZ);
//        FreeSurface_Wave_Kernel(Vector3(r,r,0),PSrc,k,G3,dX,dY,dZ);
//        qDebug() << r << G1.real() << G1.imag();// << G2.real() << G2.imag() << G3.real() << G3.imag();

        // Variation with depth

        Real z = -R*i/(NT-1);
        FreeSurface_Wave_Kernel(Vector3(2,0,z),PSrc,k,G1,dX,dY,dZ);
    }

    // Notes:
    // Imaginary part of wave kernel scales with depth as exp(-knu): Increasing k-> greater decay with depth
    // Imaginary part of wave kernel scales with radius as J0(kr):  Increasing k-> higher wavenumber.
    // ---> Select kD<<1 to ensure not significant variation of influence over panel.
    // At very low wave numbers, the J0 function is almost a constant, so that the quadrature rule produces very low errors
    // Wave number decays greatly with z, so for error checks it is better to use relative error for error metric (OR scale)
    // Qadrature is less sensitive to choice of quadrature rule, as the function varies little...
}

//--- Sub panelling (brute force num integration)

void  FlatWaveQuadPanel::CInf_SingleDoubleLayerBruteForce(const Vector3 &PGlob, CReal &S, CReal &D)
{
    // Calculate both single and double layer potentials using brute force subpanneling
    // Note ONLY TESTING: This assumes the panel normal is +z!

    int NSP = 70;
    Real DX = 2.0/NSP;

    // Subpanelling (openMP)
    CReal OS[NSP], OD[NSP];
    for (int i=0; i<NSP; i++){
        OS[i] = CReal(0,0);
        OD[i] = CReal(0,0);
    }

    OpenMPfor
    for (int i=0; i<NSP; i++){
        for (int j=0; j<NSP; j++){
            Real x = -1.0+(0.5+1.0*i)*DX;
            Real y = -1.0+(0.5+1.0*j)*DX;
            Vector3 N(0,0,1);
            Vector3 S(x,y,ZP);
            CReal G, dX, dY, dZ;
            FreeSurface_Wave_Kernel(PGlob,S,k,G,dX,dY,dZ);
//            FreeSurface_L(PGlob,S,k,G,dX,dY,dZ);         // HACK!!! ONLY L Terms!!!

            OS[i] += DX*DX*G;
            OD[i] += DX*DX*(dX*N(0)+dY*N(1)+dZ*N(2));
        }
    }

    CReal OST(0,0), ODT(0,0);
    for (int i=0; i<NSP; i++){
        OST += OS[i];
        ODT += OD[i];
    }

    // Set outputs
    S = OST;
    D = ODT;
}

}
