//-----------------------------------------------------------------------------
//-------------------------Kernel Functions------------------------------------
//-----------------------------------------------------------------------------

#include "BEMUse_Math_Types.h"

namespace BEMUse
{
//--- Source kernel

inline Real Source_Kernel(const Vector3 &S, const Vector3 &D)
{
    Vector3 R = D-S;
    Real RN = R.norm();
    if (RN!=0)  return -FourPIinv/RN;
    else        return 0.0;
}

//--- Dipole Kernel

inline Real Dipole_Kernel(const Vector3 &S, const Vector3 &D, const Vector3 &A)
{
    Vector3 r = D-S;
    Real RN = r.norm();
    if (RN!=0)  return A.dot(r)/RN/RN/RN*FourPIinv;
    else        return 0.0;
}

//--- Helper functions

static Real     BesselJ0(const Real &xx) {return std::tr1::cyl_bessel_j(0,xx);}
static Real     BesselJ1(const Real &xx) {return std::tr1::cyl_bessel_j(1,xx);}
static Real     BesselY0(const Real &xx) {return std::tr1::cyl_neumann(0,xx);}
static Real     BesselY1(const Real &xx) {return std::tr1::cyl_neumann(1,xx);}
inline Real     StruveH0(const Real &xx)
{
    // Computes the struve function of order zero

    Real SH0;

    if (xx <= 3.0)
    {
        Real yy	=	pow(xx/3.0,2);

        Real P0	=	+1.909859164;
        Real P1	=	-1.909855001;
        Real P2	=	+0.687514637;
        Real P3	=	-0.126164557;
        Real P4	=	+0.013828813;
        Real P5	=	-0.000876918;

        SH0	= P0+(P1+(P2+(P3+(P4+P5*yy)*yy)*yy)*yy)*yy;
        SH0 *= xx/3.0;
    }
    else
    {
        Real yy	=	pow(3.0/xx,2);

        Real a0	=	0.99999906;
        Real a1	=	4.77228920;
        Real a2	=	3.85542044;
        Real a3	=	0.32303607;

        Real b1	=	4.88331068;
        Real b2	=	4.28957333;
        Real b3	=	0.52120508;

        Real c1	=	2.0*(a0	+	(a1+(a2+a3*yy)*yy)*yy);
        Real c2	=	PI*xx*(1.0+	(b1+(b2+b3*yy)*yy)*yy);

        SH0	=	c1/c2 +	BesselY0(xx);
    }

    return SH0;
}
inline Real     StruveH1(const Real &xx)
{
    // Computes the derivative of the struve function of order zero

    Real SH1;

    if (xx <= 3.0)
    {
        Real yy	=	pow(xx/3.0,2);

        Real P1	=	+1.909859286;
        Real P2	=	-1.145914713;
        Real P3	=	+0.294656958;
        Real P4	=	-0.042070508;
        Real P5	=	+0.003785727;
        Real P6	=	-0.000207183;

        SH1	=	(P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy;
    }
    else
    {
        Real yy	=	pow(3.0/xx,2);

        Real a0	=	1.00000004;
        Real a1	=	3.92205313;
        Real a2	=	2.64893033;
        Real a3	=	0.27450895;

        Real b1	=	3.81095112;
        Real b2	=	2.26216956;
        Real b3	=	0.10885141;

        Real c1	=	2.0*(a0	+	(a1+(a2+a3*yy)*yy)*yy);
        Real c2	=	PI*(1.0	+	(b1+(b2+b3*yy)*yy)*yy);

        SH1	=	c1/c2	+	BesselY1(xx);
    }

    return SH1;

}

//--- Auxiliary functions: Havelock Kernel

inline Real GF_FuncA(const Real &tt)
{
    // Polynomial. Reference [2] eq. 14.a
//    Real A[10] = {1.21, -13.328, +215.896, -1763.96, +8418.94, -24314.21, +42002.57, -41592.9, 21859.0, -4838.6};
//    Real S = 0;
//    for (int i=0; i<10; i++) S += A[i]*pow(tt,i);
//    return S;

    Real S =1.21,           t = tt;
    S +=    -13.328*t;      t *= tt;
    S +=    215.896*t;      t *= tt;
    S +=    -1763.96*t;     t *= tt;
    S +=    8418.94*t;      t *= tt;
    S +=   -24314.21*t;     t *= tt;
    S +=    42002.57*t;     t *= tt;
    S +=    -41592.9*t;     t *= tt;
    S +=    21859.0*t;      t *= tt;
    S +=    -4838.6*t;      t *= tt;
    return S;
}

inline Real GF_FuncB(const Real &tt)
{
    // Polynomial. Reference [2] eq. 14.b
//    Real B[10] = {+0.938, +5.737, -67.92, +796.534, -4780.77, +17137.74, -36618.81, +44894.06, -29030.24, +7671.22};
//    Real S = 0;
//    for (int i=0; i<10; i++) S += B[i]*pow(tt,i);
//    return S;

//    Real B[10] = {+0.938, +5.737, -67.92, +796.534, -4780.77, +17137.74, -36618.81, +44894.06, -29030.24, +7671.22};
    Real S =0.938,          t = tt;
    S +=    5.737*t;        t *= tt;
    S +=    -67.92*t;       t *= tt;
    S +=    796.534*t;      t *= tt;
    S +=    -4780.77*t;     t *= tt;
    S +=    17137.74*t;     t *= tt;
    S +=    -36618.81*t;    t *= tt;
    S +=    44894.06*t;     t *= tt;
    S +=    -29030.24*t;    t *= tt;
    S +=    7671.22*t;      t *= tt;
    return S;
}

inline Real GF_FuncC(const Real &tt)
{
    // Polynomial. Reference [2] eq. 14.c
//    Real C[8] = {+1.268, -9.747, +209.653, -1397.89, +5155.67, -9844.35, +9136.4, -3272.62};
//    Real S = 0;
//    for (int i=0; i<8; i++) S += C[i]*pow(tt,i);
//    return S;

    Real S =1.268,          t = tt;
    S +=    -9.747*t;       t *= tt;
    S +=     209.653*t;     t *= tt;
    S +=    -1397.89*t;     t *= tt;
    S +=    5155.67*t;      t *= tt;
    S +=    -9844.35*t;     t *= tt;
    S +=    9136.4*t;       t *= tt;
    S +=    -3272.62*t;     t *= tt;
    return S;
}

inline Real GF_FuncD(const Real &tt)
{
    // Polynomial. Reference [2] eq. 14.c
//    Real D[10] = {+0.632,-40.97,+667.16,-6072.07,+31127.39,-96293.05,+181856.75, -205690.43,+128170.2,-33744.6};
//    Real S = 0;
//    for (int i=0; i<10; i++) S += D[i]*pow(tt,i);
//    return S;
    Real S =0.632,          t = tt;
    S +=    -40.97*t;       t *= tt;
    S +=    667.16*t;       t *= tt;
    S +=    -6072.07*t;     t *= tt;
    S +=    31127.39*t;     t *= tt;
    S +=   -96293.05*t;     t *= tt;
    S +=    181856.75*t;    t *= tt;
    S +=   -205690.43*t;    t *= tt;
    S +=    128170.2*t;     t *= tt;
    S +=   -33744.6*t;      t *= tt;
    return S;
}

inline Real GF_dFuncA(const Real &tt)
{
    // Polynomial. Reference [2] eq. 16.a
//    Real A[10] = {+2.948,-24.53,+249.69,-754.85,-1187.71,+16370.75,-48811.41,+68220.87,-46688.0,+12622.25};
//    Real S = 0;
//    for (int i=0; i<10; i++) S += A[i]*pow(tt,i);
//    return S;
    Real S =2.948,          t = tt;
    S +=    -24.53*t;       t *= tt;
    S +=    249.69*t;       t *= tt;
    S +=    -754.85*t;      t *= tt;
    S +=    -1187.71*t;     t *= tt;
    S +=    16370.75*t;     t *= tt;
    S +=    -48811.41*t;    t *= tt;
    S +=    68220.87*t;     t *= tt;
    S +=    -46688.0*t;     t *= tt;
    S +=    12622.25*t;     t *= tt;
    return S;
}

inline Real GF_dFuncB(const Real &tt)
{
    // Polynomial. Reference [2] eq. 16.b
//    Real B[10] = {+1.11,+2.894,-76.765,+1565.35,-11336.19,+44270.15,-97014.11,+118879.26,-76209.82,+19923.28};
//    Real S = 0;
//    for (int i=0; i<10; i++) S += B[i]*pow(tt,i);
//    return S;
    Real S =1.11,           t = tt;
    S +=    2.894*t;        t *= tt;
    S +=    -76.765*t;      t *= tt;
    S +=    1565.35*t;      t *= tt;
    S +=    -11336.19*t;    t *= tt;
    S +=    44270.15*t;     t *= tt;
    S +=    -97014.11*t;    t *= tt;
    S +=   118879.26*t;     t *= tt;
    S +=   -76209.82*t;     t *= tt;
    S +=    19923.28*t;     t *= tt;
    return S;
}

inline Real GF_dFuncC(const Real &tt)
{
    // Polynomial. Reference [2] eq. 16.c
//    Real C[6] = {+14.19,-148.24,+847.8,-2318.58,+3168.35,-1590.27};
//    Real S = 0;
//    for (int i=0; i<6; i++) S += C[i]*pow(tt,i);
//    return S;
    Real S =14.19,          t = tt;
    S +=    -148.24*t;      t *= tt;
    S +=    847.8*t;        t *= tt;
    S +=    -2318.58*t;     t *= tt;
    S +=    3168.35*t;      t *= tt;
    S +=   -1590.27*t;      t *= tt;
    return S;
}

inline Real GF_Func_L0(Real *GF_Vars, const Real &hh, const Real &vv)
{
    // Flow part of GF Kernel
    Real PP = log(0.5*(GF_Vars[0]-vv)) + EUL - 2.0*pow(GF_Vars[0],2);
    PP *= exp(vv);
    PP += pow(GF_Vars[0],2) - vv;

    // Replace: GF_Func_Lp function

    Real A = GF_FuncA(GF_Vars[4]);
    Real B = GF_FuncB(GF_Vars[4]);
    Real C = GF_FuncC(GF_Vars[4]);
    Real D = GF_FuncD(GF_Vars[4]);

    Real RR = (1.0-GF_Vars[2])*A;
    RR -= GF_Vars[2]*B;
    RR -= GF_Vars[1]*C/(1.0+6.0*GF_Vars[1]*GF_Vars[4]*(1.0-GF_Vars[4]));
    RR += GF_Vars[2]*(1.0-GF_Vars[2])*D;

    Real GF_Func_Lp = GF_Vars[4]*pow(1.0-GF_Vars[4],3)*RR;

    Real GF_Func_L0 = 2.0*PP/(1.0+pow(GF_Vars[0],3)) + 2.0*GF_Func_Lp;

    return GF_Func_L0;
}

inline CReal GF_Func_W(Real *GF_Vars, const Real &hh, const Real &vv)
{
    // Wave part of GF Kernel

    Real H0 = StruveH0(hh);
    Real J0 = BesselJ0(hh);

    Real E = 2.0*PI*exp(vv);
    CReal GFW_Out = CReal(E*H0,-E*J0);

    return GFW_Out;
}

inline Real GF_Func_Ls(Real *GF_Vars, const Real &hh, const Real &vv)
{
    // This calculates the spatial derivates of the Green. This corresponds to the dipole component
    // of the Greens function wave component

    Real PS = (GF_Vars[2]+hh)/(GF_Vars[0]-vv);
    PS -= 2.0*GF_Vars[2];
    PS += 2.0*GF_Vars[0]*exp(vv);
    PS -= hh;

    Real QS = exp(-GF_Vars[0])*(1.0-GF_Vars[2]);
    QS *= 1.0 + GF_Vars[0]/(1.0+pow(GF_Vars[0],3));

    // Replace: GF_Func_Lsp function

    Real A = GF_dFuncA(GF_Vars[4]);
    Real B = GF_dFuncB(GF_Vars[4]);
    Real C = GF_dFuncC(GF_Vars[4]);

    Real RR	= GF_Vars[2]*A;
    RR	-=	(1.0-GF_Vars[1])*B;
    RR	+=	GF_Vars[2]*(1.0-GF_Vars[2])*GF_Vars[4]*(1.0-2.0*GF_Vars[4])*C;

    Real Lsp =	GF_Vars[4]*pow(1.0-GF_Vars[4],3)*RR;

    Real GF_Func_Ls = 2.0*PS/(1.0+pow(GF_Vars[0],3)) - 4.0*QS + 2.0*Lsp;

    return GF_Func_Ls;
}

inline CReal GF_Func_Wh(Real *GF_Vars, const Real &hh, const Real &vv)
{
    // Wave part of GF Kernel

    Real H1 = StruveH1(hh);
    Real J1 = BesselJ1(hh);

    Real E = 2.0*PI*exp(vv);
    CReal GFWh_Out = CReal(E*(2.0/PI-H1),E*J1);

    return GFWh_Out;
}

//--- Liang/Noblesse implementation of Havelock kernel

inline void HavelockGF(const Real &dx, const Real &dy, const Real &vv, CReal &GF, CReal &DGFX, CReal &DGFY, CReal &DGFZ)
{
    //  This calculates the wave component of the greens functions for the case of infinite depth

    //    Code converted to C++ from https://www.researchgate.net/profile/Francis_Noblesse. Thanks guys!

    //    This subroutine evaluates the free-surface term
    //    of Green function and its spatial derivatives
    //    for the wave radiation-diffraction problem.
    //    !
    //    The local-flow component is approximated by
    //    mean of the global approximations [1].
    //    The computations reported in [2] provides
    //    strong evidence that the global approximations
    //    are sufficiently accurate to compute linear
    //    and second-order wave loads in practice.
    //    !
    //    It should be noted that the Rankine source
    //    term -1/d appeared in L_z given by (8a) in
    //    Ref. [2] is not evaluated here. These Rankine
    //    source components are required to evaluate
    //    in another routine.

    //    For any questions, please contact: lianghuistar@gmail.com

    //    [1]	H. Wu, C. Zhang, Y. Zhu, W. Li, D. Wan, F. Noblesse,
    //        A global approximation to the Green function for
    //        diffraction radiation of water waves,
    //        Eur. J. Mech. B Fluids 65 (2017) 54-64.

    //    [2]	H. Liang, H. Wu, F. Noblesse,
    //        Validation of a global approximation for
    //        wave diffraction-radiation in deep water,
    //        Appl. Ocean Res. 74 (2018) 80-86.

    Real hh = sqrt(pow(dx,2) + pow(dy,2));

    // Replace: GF_DivParameters function

    Real GF_Vars[5];

    GF_Vars[0] = sqrt(pow(hh,2) + pow(vv,2));   // GF_dd
    GF_Vars[1] = -vv/GF_Vars[0];                // GF_alpha
    GF_Vars[2] = hh/GF_Vars[0];                 // GF_beta
    GF_Vars[3] = hh/(GF_Vars[0]-vv);            // GF_sigma
    GF_Vars[4] = GF_Vars[0]/(1.0+GF_Vars[0]);   // GF_rho

    // First calculate the greens function
    CReal GF_L_Part = CReal(GF_Func_L0(GF_Vars,hh,vv),0.0);
    CReal GF_W_Part = GF_Func_W(GF_Vars,hh,vv);

    GF = GF_L_Part + GF_W_Part;

    // Now calculate spatial derivs
    CReal GFh_L = CReal(GF_Func_Ls(GF_Vars,hh,vv),0.0);
    CReal GFh_W = GF_Func_Wh(GF_Vars,hh,vv);
    CReal GFh = GFh_L + GFh_W;

    if (hh > 1.0e-6)
    {
        DGFX = GFh*dx/hh;
        DGFY = GFh*dy/hh;
    }
    else
    {
        DGFX = CReal(0.0,0.0);
        DGFY = CReal(0.0,0.0);
    }

    DGFZ = GF;
}

}
