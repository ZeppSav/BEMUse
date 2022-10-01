/****************************************************************************
    Airfoil Profiles Class
    Copyright (C) 2019 Joseph Saverin j.saverin@tu-berlin.de
    Naca profiles from Ludwig Koennecke

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

    Information on class:

    -> Simply creates the coordinates of the chosen airfoils. All members are static

*****************************************************************************/

#ifndef AIRFOIL_PROFILES_H
#define AIRFOIL_PROFILES_H

#include "../BEMUse_Math_Types.h"

namespace BEMUse
{

inline void Karman_Trefftz(CReal &Centre, Real &TE_Angle, int Disc, StdVector &X, StdVector &Y)
{
    // Creates a 2D Profile of a Karman-Trefftz airfoil
    // This is a copy of the function described in the document:
    // Kerwin_2001: "13.04 Lecture Notes Hydrofoils and Propellors"

    // Calc the "appropriate" circle radius for the given centre
    Real Lambda = 2.0 - TE_Angle/180.0;
    Real RCSQ = pow(1.0-Centre.real(),2) + pow(Centre.imag(),2);
    Real RC = sqrt(RCSQ);

    // Create Z and Zeta plane coords
    for (int i=0; i<Disc; i++)
    {
        Real Theta = i*2.0*PI/Disc;
        Real XP = Centre.real() + RC*cos(Theta);
        Real YP = Centre.imag() + RC*sin(Theta);

        // Map positions to ZETA Plane
        CReal ZP11 (XP + 1.0, YP);
        CReal ZM11 (XP - 1.0, YP);
        CReal ZP1 =  pow(ZP11,Lambda);
        CReal ZM1 =  pow(ZM11,Lambda);
        CReal Zeta = Lambda*(ZP1+ZM1)/(ZP1-ZM1);
        X.push_back(Zeta.real());
        Y.push_back(Zeta.imag());
    }
    X.push_back(X[0]);
    Y.push_back(Y[0]);
}

inline void Naca_4digit(int Disc, int N1, int N2, int N3, Vector &X, Vector &Y)
{

    // Set the input parameters here

    Real A=N1;               //First Number of Naca Profile - Profile camber in % of the length
    Real B=N2;               //Second Number of Naca Profile - Position of the camber in 1/10 of the length
    Real CD=N3;             //Last two Numbers of Naca Profile - thickness in % of the length

    Real in_AoA = 0;       // Angle of Attack in °
//    Real AoA = -in_AoA*(PI/180.0);     //in Radian

    Real M=(A/100.0);    //maximum camber in %
    Real P=(B/10.0);     // position of the max camber in tenth
    Real T=(CD/100.0);   //thickness in % of the chord

    Real beta;
    Real x;
    Real yc;
    Real dyc_dx;
    Real yt;
    Real phi;

    // Coefficients for the calculation of the thickness thats added to the camber
    Real a0=0.2969;
    Real a1=-0.126;
    Real a2=-0.3516;
    Real a3=0.2843;
    Real a4=-0.1036; // -0.1015 for a finite thick T.E.

//Loop to calculate the Grid Points of the Naca Profile

    for (int i=0; i<(Disc*0.5); i++)
    {
        //cosinus distribution
//        beta=pow((Disc-1),-1)*2*i*PI;       //to get a unifrom distance between the x values
//        x=(1+cos(beta))/2;             //x values on the chord

          //linear and cosinus distribution
        beta=pow((Disc-1),-1)*i*2*PI;       //to get a uniform distance between the x values
        if (i<=(Disc*0.25)){
            x=1-(i/(Disc*0.5));
        }
        else{
            x=(1+cos(beta))/2;
        }


        //Camberline before max. camber (0 <= x < P)
        if(0<=x && x<P)
        {
            yc = (M/(pow(P,2)))*(2*P*x-pow(x,2));
            dyc_dx=((2*M)/(pow(P,2)))*(P-x);
        }
        //Camberline after max.camber (P <= x <=1)
        else if (P<=x && x<= 1)
        {
           yc=(M/(pow((1-P),2)))*(1-2*P+2*P*x-pow(x,2));
           dyc_dx=((2*M)/(pow(1-P,2)))*(P-x);
        }
        else {
            yc=0;
            dyc_dx=0;
        }

        yt=(T/0.2)*(a0*pow(x,0.5)+a1*x+a2*pow(x,2)+a3*pow(x,3)+a4*pow(x,4));        //calculates the thickness

        phi= atan(dyc_dx);

        X(i)=x-yt*sin(phi);
        Y(i)=yc+yt*cos(phi);

        X((Disc-1)-i)=x+yt*sin(phi);
        Y((Disc-1)-i)=yc-yt*cos(phi);
    }
    Y(0)=0;
    Y(Disc-1)=0;
}

inline void Naca_5digit(int Disc, int N1, int N2, int N3, int N4, Vector &X, Vector &Y)
{

    // Set the input parameters here

    Real A=N1;               //First Number of Naca Profile - deisgn coefficient of lift * 0.15 (r,k1,k2k1 are for Cl of 0.3 in this program)
    Real B=N2;               //Second Number of Naca Profile - Position of max. camber in 1/20 of the chord
    Real C=N3;               //Third Numbers of Naca Profile - 0 = normal camber line, 1 reflex camber line
    Real DE=N4;             //Last two numbers of Nacka Profile - thickness in % of the chord

    Real in_AoA = 0;       // Angle of Attack in °
    Real AoA = in_AoA*(PI/180);     //in Radian

    Real P=B/20;        //pos. of max camber of chord
    Real T=(DE/100);   //thickness of the chord

    Real beta;
    Real x;
    double yc;
    double dyc_dx;
    double yt;
    double phi;

    Real r;
    Real k1;
    Real k2k1;

    // Coefficients for the calculation of the thickness thats added to the camber
    double a0=0.2969;
    double a1=-0.126;
    double a2=-0.3516;
    double a3=0.2843;
    double a4=-0.1015;

    if (C==0)
    {      //normal camber line
        for (int i=0; i<(0.5*Disc); i++)
        {
            if(B==1){
                r=0.0580;
                k1=361.400;
            }
            else if(B==2){
                r=0.1260;
                k1=51.640;
            }
            else if(B==3){
                r=0.2025;
                k1=15.957;
            }
            else if(B==4){
                r=0.2900;
                k1=6.643;
            }
            else if(P==5){
                r=0.3910;
                k1=3.230;
            }
            else{
                r=0;
                k1=0;
            }

            beta=pow((Disc-1),-1)*2*i*PI;       //to get a unifrom distance between the x values
            x=((1-cos(beta))/2);                //x values on the chord

            //Camberline before max. camber (0 <= x < P)
            if(0<=x && x<P)
            {
                yc = (k1/6)*(pow(x,3)-3*r*pow(x,2)+pow(r,2)*(3-r)*x);
                dyc_dx=(k1/6)*(3*pow(x,2)-6*r*x+pow(r,2)*(3-r));
            }
            //Camberline after max.camber (P <= x <=1)
            else if (P<=x && x<= 1)
            {
               yc=((k1*pow(r,3))/6)*(1-x);
               dyc_dx=-((k1*pow(r,3))/6);
            }
            else {
                yc=0;
                dyc_dx=0;
            }

            yt=(T/0.2)*(a0*pow(x,0.5)+a1*x+a2*pow(x,2)+a3*pow(x,3)+a4*pow(x,4));        //calculates the thickness

            phi= atan(dyc_dx);

            Real Xoo=x-yt*sin(phi);
            Real Yoo=yc+yt*cos(phi);

            Real Xou=x+yt*sin(phi);
            Real You=yc-yt*cos(phi);

            //For Roation of the Profile

            X(i) = cos(AoA)*Xoo - sin(AoA)*Yoo;
            Y(i) = sin(AoA)*Xoo + cos(AoA)*Yoo;

            X((Disc-1)-i) = cos(AoA)*Xou - sin(AoA)*You;
            Y((Disc-1)-i) = sin(AoA)*Xou + cos(AoA)*You;
        }
    }
    else if (C==1)
    {      //reflex camber line
        for (int i=0; i<(0.5*Disc); i++)
        {
            if(B==2){
                r=0.1300;
                k1=51.990;
                k2k1=0.000764;
            }
            else if(B==3){
                r=0.2170;
                k1=15.793;
                k2k1=0.00677;
            }
            else if(B==4){
                r=0.3180;
                k1=6.520;
                k2k1=0.0303;
            }
            else if(B==5){
                r=0.4410;
                k1=3.191;
                k2k1=0.1355;
            }
            else{
                r=0;
                k1=0;
                k2k1=0;
            }

            beta=pow((Disc-1),-1)*2*i*PI;       //to get a unifrom distance between the x values
            x=((1-cos(beta))/2);                //x values on the chord

            //Camberline before max. camber (0 <= x < P)
            if(0<=x && x<P)
            {
                yc=(k1/6)*(pow((x-r),3)-k2k1*pow((1-r),3)*x-pow(r,3)*x+(pow(r,3)));
                dyc_dx=(k1/6)*(3*pow((x-r),2)-k2k1*pow((1-r),3)-pow(r,3));
            }
            //Camberline after max.camber (P <= x <=1)
            else if (P<=x && x<= 1)
            {
               yc=(k1/6)*(k2k1*pow((x-r),3)-k2k1*pow((1-r),3)*x-pow(r,3)*x+(pow(r,3)));
               dyc_dx=(k1/6)*(3*k2k1*pow((x-r),2)-k2k1*pow((1-r),3)-pow(r,3));
            }
            else {
                yc=0;
                dyc_dx=0;
            }

            yt=(T/0.2)*(a0*pow(x,0.5)+a1*x+a2*pow(x,2)+a3*pow(x,3)+a4*pow(x,4));        //calculates the thickness

            phi= atan(dyc_dx);

            Real Xoo=x-yt*sin(phi);
            Real Yoo=yc+yt*cos(phi);

            Real Xou=x+yt*sin(phi);
            Real You=yc-yt*cos(phi);

            //For Roation of the Profile

            X(i) = cos(AoA)*Xoo - sin(AoA)*Yoo;
            Y(i) = sin(AoA)*Xoo + cos(AoA)*Yoo;

            X((Disc-1)-i) = cos(AoA)*Xou - sin(AoA)*You;
            Y((Disc-1)-i) = sin(AoA)*Xou + cos(AoA)*You;
        }
    }
    else std::cout << "Wrong input of the third Naca number. It's value should be 0 or 1.\n";

}

inline void Normalize(int Disc, StdVector &X, StdVector &Y)
{
    // This simply scales the airfoil to lie between 0 <= X <= 1

    // Now Normalize
//    Matrix_ID Min_C, Min_R, Max_C, Max_R;
    Real X_min  = *std::min_element(X.begin(),X.end());
    Real X_max  = *std::max_element(X.begin(),X.end());
//    Real X_min = X.minCoeff(&Min_C,&Min_R);
//    Real X_max = X.maxCoeff(&Max_C,&Max_R);
    Real Scale = 1.0/(X_max-X_min);

    // Shift X
    for (int C=0; C<Disc; C++)  X[C] -= X_min;
    // Scale
    for (int C=0; C<Disc; C++)
    {
        X[C] *= Scale;
        Y[C] *= Scale;
    }
}

inline void Seal_Sharp_TE(int Disc, Vector &X, Vector &Y)
{
    // This is simply a function to check if the TE is sharp.
    // If not then an additional node is added

    // NB: This is only called if we are SURE or DICTATE that the TE is sharp.
    // Furthermore: This ONLY works on the KT profile!

    Vector3 P1; P1 << X(0), Y(0), 0.0;
    Vector3 P2; P2 << X(Disc-1), Y(Disc-1), 0.0;

    Vector3 RelPos = P1-P2;

    if (RelPos.norm() > 1e-6)
    {
        std::cout << "Seal_Sharp_TE: Sealing Trailing edge.\n";

        Vector X2 = Vector::Zero(Disc+1);
        Vector Y2 = Vector::Zero(Disc+1);

        X2.block(1,0,Disc,1) = X;
        Y2.block(1,0,Disc,1) = Y;

        X2(0) = X(Disc-1);
        Y2(0) = Y(Disc-1);

        Disc++;

        X = X2;
        Y = Y2;
    }
}

}

#endif // AIRFOIL_PROFILES_H
