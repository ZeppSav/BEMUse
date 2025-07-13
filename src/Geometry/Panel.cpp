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

//--- Linear strength distribution

void FlatSourceTriPanel::Inf_SingleLayer_LinDist(const Vector3 &P, Real &S1, Real &S2, Real &S3)
{
    // This implementation is from Letournel (doi: )

    // Calculate geometric parameters.
    Vector3 CLoc = Centroid()->Position_Global();
    Matrix3 OLoc = Centroid()->OrientMat_Global();

    Vector3 P0 = OLoc * (Geo->Nodes[0]->Position_Global() - CLoc);
    Vector3 P1 = OLoc * (Geo->Nodes[1]->Position_Global() - CLoc);
    Vector3 P2 = OLoc * (Geo->Nodes[2]->Position_Global() - CLoc);
    Vector3 PG = OLoc * (Geo->Centroid->Position_Global() - CLoc);

    // Normal of source panel
    Vector3 P_Normal = Centroid()->Z_Axis_Global();

    // Receiver node
    Vector3 PL = OLoc * (P - CLoc);

    Matrix3 AM;
    Matrix3 AB;

    AM.row(0) = PL - P0;
    AM.row(1) = PL - P1;
    AM.row(2) = PL - P2;

    // Calculate edge vectors of the triangle
    AB.row(0) = P1 - P0;
    AB.row(1) = P2 - P1;
    AB.row(2) = P0 - P2;

    // -----------------------
    // Delta Calculation (see Letournel)
    // -----------------------
    Vector3 P0P1 = AB.row(0);
    Vector3 P0P2 = -AB.row(2);

    Real P0P1_sq = P0P1.squaredNorm();
    Real P0P2_sq = P0P2.squaredNorm();

    Real P0P1_P0P2_dot = P0P1.dot(P0P2);

    Real Delta = P0P1_sq * P0P2_sq - P0P1_P0P2_dot * P0P1_P0P2_dot;

    // -----------------------
    // A B Calculation (see Letournel)
    // -----------------------
    Vector3 A = P0P2_sq * P0P1 - P0P1_P0P2_dot * P0P2;
    Vector3 B = -P0P1_P0P2_dot * P0P1 + P0P1_sq * P0P2;

    // -----------------------
    // Sigma Calculation (see Letournel)
    // -----------------------
    Matrix3 Sigma;
    Sigma.row(0) = -(A + B);
    Sigma.row(1) = A;
    Sigma.row(2) = B;
    Sigma *= 1. / Delta;

    // -----------------------
    // S_sigma Calculation (see Letournel)
    // -----------------------
    Vector3 R_k = AM.rowwise().norm();
    Vector3 d_k;
    Vector3 N_k;
    Vector3 N_k_1;
    Vector3 D_k_1;
    Vector3 D_k;

    Real S_sigma = 0.;

    Real Z = (PL - PG).dot(P_Normal);
    Real Z_abs = abs(Z);

    // -----------------------
    // Line Integral I_sigma
    // -----------------------
    Vector3 I_sigma(0., 0., 0.);
    Vector3 a_lineIntegral;
    Vector3 b_lineIntegral;
    Vector3 q0;
    Vector3 q1;
    Vector3 K2;

    Real Temp;

    for (int j = 0; j < 3; ++j) {
        Real AM_sq_norm = AM.row(j).squaredNorm();
        Real AB_AM_dot = AB.row(j).dot(AM.row(j));
        Real AB_sq_norm = AB.row(j).squaredNorm();
        Real AB_norm = sqrt(AB_sq_norm);

        // S_sigma
        int next = (j + 1) % 3;

        N_k(j) = 2. * AM.row(j).dot(P_Normal.cross(AB.row(j)));
        d_k(j) = AB_norm;

        N_k_1(j) = R_k(next) + R_k(j) + d_k(j);
        D_k_1(j) = R_k(next) + R_k(j) - d_k(j);
        D_k(j) = (R_k(next) + R_k(j)) * (R_k(next) + R_k(j)) - d_k(j) * d_k(j) + 2 * Z_abs * (R_k(next) + R_k(j));

        S_sigma += 0.5 * N_k(j) / d_k(j) * log(N_k_1(j) / D_k_1(j)) - 2 * Z_abs * atan(N_k(j) / D_k(j));

        // I_sigma

        K2(j) = AM_sq_norm - AB_AM_dot * AB_AM_dot / AB_sq_norm;

        Real K2_sqrt = sqrt(K2(j));

        q0(j) = -AB_AM_dot / (AB_norm * K2_sqrt);
        q1(j) = (AB_sq_norm - AB_AM_dot) / (AB_norm * K2_sqrt);

        a_lineIntegral(j) = asinh(q0(j));
        b_lineIntegral(j) = asinh(q1(j));

        Real denom = 2.0 * AB_norm;
        Real sinh_2a = sinh(2.0 * a_lineIntegral(j));
        Real sinh_2b = sinh(2.0 * b_lineIntegral(j));

        Temp = K2(j) / denom * (b_lineIntegral(j) - a_lineIntegral(j) + 0.5 * (sinh_2b - sinh_2a));

        I_sigma += Temp * (P_Normal.cross(AB.row(j)));
    }

    // -----------------------
    // SRC and DIP Influence
    // -----------------------
    Vector3 I_single;
    Vector3 P_I(1., 1., 1.);
    P_I = P_I / 3.;

    I_single = (P_I + Sigma * (PL - PG)) * S_sigma - Sigma * I_sigma;

    S1 = I_single(0) * FourPIinv;
    S2 = I_single(1) * FourPIinv;
    S3 = I_single(2) * FourPIinv;
}

void FlatSourceTriPanel::Inf_DoubleLayer_LinDist(const Vector3 &P, Real &D1, Real &D2, Real &D3)
{
    // This implementation is from Letournel (doi: )

    // Calculate geometric parameters.
    Vector3 CLoc = Centroid()->Position_Global();
    Matrix3 OLoc = Centroid()->OrientMat_Global();

    // std::cout << "CLoc " << CLoc << " complete.\n" << std::endl;
    // std::cout << "OLoc " << OLoc << " complete.\n" << std::endl;

    Vector3 P0 = OLoc * (Geo->Nodes[0]->Position_Global() - CLoc);
    Vector3 P1 = OLoc * (Geo->Nodes[1]->Position_Global() - CLoc);
    Vector3 P2 = OLoc * (Geo->Nodes[2]->Position_Global() - CLoc);
    Vector3 PG = OLoc * (Geo->Centroid->Position_Global() - CLoc);

    // Normal of source panel
    Vector3 P_Normal = Centroid()->Z_Axis_Global();

    // Receiver node
    Vector3 PL = OLoc * (P - CLoc);

    Matrix3 AM;
    Matrix3 AB;

    AM.row(0) = PL - P0;
    AM.row(1) = PL - P1;
    AM.row(2) = PL - P2;

    // Calculate edge vectors of the triangle
    AB.row(0) = P1 - P0;
    AB.row(1) = P2 - P1;
    AB.row(2) = P0 - P2;

    // -----------------------
    // Delta Calculation (see Letournel)
    // -----------------------
    Vector3 P0P1 = AB.row(0);
    Vector3 P0P2 = -AB.row(2);

    Real P0P1_sq = P0P1.squaredNorm();
    Real P0P2_sq = P0P2.squaredNorm();

    Real P0P1_P0P2_dot = P0P1.dot(P0P2);

    Real Delta = P0P1_sq * P0P2_sq - P0P1_P0P2_dot * P0P1_P0P2_dot;

    // -----------------------
    // A B Calculation (see Letournel)
    // -----------------------
    Vector3 A = P0P2_sq * P0P1 - P0P1_P0P2_dot * P0P2;
    Vector3 B = -P0P1_P0P2_dot * P0P1 + P0P1_sq * P0P2;

    // -----------------------
    // Sigma Calculation (see Letournel)
    // -----------------------
    Matrix3 Sigma;
    Sigma.row(0) = -(A + B);
    Sigma.row(1) = A;
    Sigma.row(2) = B;

    Sigma *= 1. / Delta;

    // -----------------------
    // S_sigma Calculation (see Letournel)
    // -----------------------
    Vector3 R_k = AM.rowwise().norm();
    Vector3 d_k;
    Vector3 N_k;
    Vector3 N_k_1;
    Vector3 D_k_1;
    Vector3 D_k;

    Real Z = (PL - PG).dot(P_Normal);
    Real Z_abs = abs(Z);

    // -----------------------
    // Line Integral I_sigma
    // -----------------------
    Vector3 q0;
    Vector3 q1;
    Vector3 K2;

    // -----------------------
    // S_mu Calculation (see Letournel)
    // -----------------------
    Real S_mu = 0.;
    Real Z_Sign;

    if (Z > 0) {
        Z_Sign = 1;
    } else {
        Z_Sign = -1;
    }

    // ------------------------
    // I_mu Calculation (see Letournel)
    // ------------------------
    Vector3 I_mu(0., 0., 0.);
    Real logterm_q;

    for (int j = 0; j < 3; ++j) {
        Real AM_sq_norm = AM.row(j).squaredNorm();
        Real AB_AM_dot = AB.row(j).dot(AM.row(j));
        Real AB_sq_norm = AB.row(j).squaredNorm();
        Real AB_norm = sqrt(AB_sq_norm);

        int next = (j + 1) % 3;

        N_k(j) = 2. * AM.row(j).dot(P_Normal.cross(AB.row(j)));
        d_k(j) = AB_norm;

        N_k_1(j) = R_k(next) + R_k(j) + d_k(j);
        D_k_1(j) = R_k(next) + R_k(j) - d_k(j);
        D_k(j) = (R_k(next) + R_k(j)) * (R_k(next) + R_k(j)) - d_k(j) * d_k(j) + 2 * Z_abs * (R_k(next) + R_k(j));

        K2(j) = AM_sq_norm - AB_AM_dot * AB_AM_dot / AB_sq_norm;

        Real K2_sqrt = sqrt(K2(j));

        q0(j) = -AB_AM_dot / (AB_norm * K2_sqrt);
        q1(j) = (AB_sq_norm - AB_AM_dot) / (AB_norm * K2_sqrt);

        // S_mu
        S_mu += 2 * Z_Sign * atan(N_k(j) / D_k(j));

        // I_mu
        Real sqrt_1_plus_q0_sq = sqrt(1. + q0(j) * q0(j));
        Real sqrt_1_plus_q1_sq = sqrt(1. + q1(j) * q1(j));

        logterm_q = log((q1(j) + sqrt_1_plus_q1_sq) / (q0(j) + sqrt_1_plus_q0_sq));

        I_mu += (AM.row(j).cross(AB.row(j))) / AB_norm * logterm_q;
    }

    Vector3 I_double;
    Vector3 P_I(1., 1., 1.);
    P_I = P_I / 3.;

    I_double = -1. * (P_I + Sigma * (PL - PG)) * S_mu - Sigma * I_mu;

    D1 = I_double(0) * FourPIinv;
    D2 = I_double(1) * FourPIinv;
    D3 = I_double(2) * FourPIinv;
}

void FlatSourceTriPanel::Inf_SingleDoubleLayer_LinDist(const Vector3 &P, Real &S1, Real &S2, Real &S3, Real &D1, Real &D2, Real &D3)
{
    // Newman, J.N., "Distributions of sources and normal dipoles over a quadrilateral panel",
    // 1986, https://doi.org/10.1007/BF00042771

    Real R[3], S[3], m[3], cs[3], sn[3], B[3], e[3], h[3], Q[3], u[3], U[3], Pn[3];
    Vector3 r[3], s[3];

    Vector3 CLoc = Centroid()->Position_Global();
    Matrix3 OLoc = Centroid()->OrientMat_Global();

    Vector3 P0 = OLoc * (Geo->Nodes[0]->Position_Global() - CLoc);
    Vector3 P1 = OLoc * (Geo->Nodes[1]->Position_Global() - CLoc);
    Vector3 P2 = OLoc * (Geo->Nodes[2]->Position_Global() - CLoc);
    Vector3 PL = OLoc * (P - CLoc);

    // if (PL(2)==0) PL(2) -= 1.e-16;          // Avoid issues with evaluation points exactly at zero.
    // Real dFac = 1.e-8;
    // if ((PL-P0).norm() < dFac || (PL-P1).norm() < dFac|| (PL-P2).norm() < dFac){
    //     S1 = 0.; S2 = 0.; S3 = 0.;
    //     D1 = 0.; D2 = 0.; D3 = 0.;
    //     return;
    // }

    r[0] = PL - P0; R[0] = r[0].norm();
    r[1] = PL - P1; R[1] = r[1].norm();
    r[2] = PL - P2; R[2] = r[2].norm();

    s[0] = P1 - P0; S[0] = s[0].norm();
    s[1] = P2 - P1; S[1] = s[1].norm();
    s[2] = P0 - P2; S[2] = s[2].norm();

    m[0] = (P1(1) - P0(1)) / (P1(0) - P0(0));
    m[1] = (P2(1) - P1(1)) / (P2(0) - P1(0));
    m[2] = (P0(1) - P2(1)) / (P0(0) - P2(0));

    //cos and sin of panel
    cs[0] = (P1(0) - P0(0)) / S[0];
    cs[1] = (P2(0) - P1(0)) / S[1];
    cs[2] = (P0(0) - P2(0)) / S[2];

    sn[0] = (P1(1) - P0(1)) / S[0];
    sn[1] = (P2(1) - P1(1)) / S[1];
    sn[2] = (P0(1) - P2(1)) / S[2];

    B[0] = (PL(0) - P0(0)) * sn[0] - (PL(1) - P0(1)) * cs[0];
    B[1] = (PL(0) - P1(0)) * sn[1] - (PL(1) - P1(1)) * cs[1];
    B[2] = (PL(0) - P2(0)) * sn[2] - (PL(1) - P2(1)) * cs[2];

    e[0] = (PL(0) - P0(0)) * (PL(0) - P0(0)) + PL(2) * PL(2);
    e[1] = (PL(0) - P1(0)) * (PL(0) - P1(0)) + PL(2) * PL(2);
    e[2] = (PL(0) - P2(0)) * (PL(0) - P2(0)) + PL(2) * PL(2);

    h[0] = (PL(0) - P0(0)) * (PL(1) - P0(1));
    h[1] = (PL(0) - P1(0)) * (PL(1) - P1(1));
    h[2] = (PL(0) - P2(0)) * (PL(1) - P2(1));

    Q[0] = log((R[0] + R[1] + S[0]) / (R[0] + R[1] - S[0]));
    Q[1] = log((R[1] + R[2] + S[1]) / (R[1] + R[2] - S[1]));
    Q[2] = log((R[2] + R[0] + S[2]) / (R[2] + R[0] - S[2]));

    U[0] = r[1].dot(s[0] / S[0]);
    U[1] = r[2].dot(s[1] / S[1]);
    U[2] = r[0].dot(s[2] / S[2]);

    u[0] = r[0].dot(s[0] / S[0]);
    u[1] = r[1].dot(s[1] / S[1]);
    u[2] = r[2].dot(s[2] / S[2]);

    Pn[0] = 0.5 * (u[0] * R[0] - U[0] * R[1] + (R[0] * R[0] - u[0] * u[0]) * Q[0]);
    Pn[1] = 0.5 * (u[1] * R[1] - U[1] * R[2] + (R[1] * R[1] - u[1] * u[1]) * Q[1]);
    Pn[2] = 0.5 * (u[2] * R[2] - U[2] * R[0] + (R[2] * R[2] - u[2] * u[2]) * Q[2]);

    Real PhiDipConst =
        -atan((m[0] * e[0] - h[0]) / (PL(2) * R[0])) + atan((m[0] * e[1] - h[1]) / (PL(2) * R[1]))
        -atan((m[1] * e[1] - h[1]) / (PL(2) * R[1])) + atan((m[1] * e[2] - h[2]) / (PL(2) * R[2]))
        -atan((m[2] * e[2] - h[2]) / (PL(2) * R[2])) + atan((m[2] * e[0] - h[0]) / (PL(2) * R[0]));

    Real PhiSrcConst = -(B[0] * Q[0] + B[1] * Q[1] + B[2] * Q[2]) - PL(2) * PhiDipConst;

    Real SrcsumX = Pn[0] * sn[0] + Pn[1] * sn[1] + Pn[2] * sn[2];
    Real SrcsumY = Pn[0] * cs[0] + Pn[1] * cs[1] + Pn[2] * cs[2];

    Real DipsumX = Q[0] * sn[0] + Q[1] * sn[1] + Q[2] * sn[2];
    Real DipsumY = Q[0] * cs[0] + Q[1] * cs[1] + Q[2] * cs[2];

    // Final constants and results
    Real DipConst = FourPIinv * PhiDipConst;
    Real SrcConst = FourPIinv * PhiSrcConst;

    Real DipX = -FourPIinv * (-PL(0) * PhiDipConst + PL(2) * DipsumX);
    Real DipY = -FourPIinv * (-PL(1) * PhiDipConst - PL(2) * DipsumY);
    Real SrcX = -FourPIinv * (-PL(0) * PhiSrcConst - SrcsumX);
    Real SrcY = -FourPIinv * (-PL(1) * PhiSrcConst + SrcsumY);

    Matrix3 Surface_Dist = Matrix3::Ones();
    Surface_Dist(0,1) = P0(0);    Surface_Dist(0,2) = P0(1);
    Surface_Dist(1,1) = P1(0);    Surface_Dist(1,2) = P1(1);
    Surface_Dist(2,1) = P2(0);    Surface_Dist(2,2) = P2(1);
    Matrix3 Surface_Interp = Surface_Dist.inverse();

    Vector3 sdP1 = Surface_Interp*Vector3(1,0,0);   // Surface distribution for case that node 1 has unit strength
    Vector3 sdP2 = Surface_Interp*Vector3(0,1,0);   // Surface distribution for case that node 2 has unit strength
    Vector3 sdP3 = Surface_Interp*Vector3(0,0,1);   // Surface distribution for case that node 3 has unit strength

    Vector3 SVec(SrcConst,SrcX,SrcY);
    Vector3 DVec(DipConst,DipX,DipY);

    S1 = sdP1.dot(SVec);
    S2 = sdP2.dot(SVec);
    S3 = sdP3.dot(SVec);
    D1 = sdP1.dot(DVec);
    D2 = sdP2.dot(DVec);
    D3 = sdP3.dot(DVec);
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
