/****************************************************************************
    CoordSys class
    Copyright (C) 2021 Joseph Saverin j.saverin@tu-berlin.de

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

    -> This defines the coordinate system type.

*****************************************************************************/

#ifndef COORDSYS_H
#define COORDSYS_H

#include "../BEMUse_Math_Types.h"

namespace BEMUse
{

static Vector3 UnitX = Vector3::UnitX();
static Vector3 UnitY = Vector3::UnitY();
static Vector3 UnitZ = Vector3::UnitZ();

class CoordSys
{
private:

    Vector3 R = Vector3::Zero();    // Position vector within the reference system
    Vector3 U = Vector3::Zero();    // Velocity vector within the reference system
    Vector3 W = Vector3::Zero();    // Angular velocity vector within the reference system
    Quat O = Quat::Identity();      // Coordinate system orientation
    CoordSys *Ref_CoordSys = nullptr;// The coordinate system within which the position and rotations are defined.

public:

    CoordSys()                      {}
    CoordSys(CoordSys *tRCS)        {Ref_CoordSys = tRCS;}
    CoordSys(CoordSys *tRCS, const Vector3 &tR) : CoordSys(tRCS)    {R = tR;}
    CoordSys(CoordSys *tRCS, const Quat &tQ, const Vector3 &tR)  : CoordSys(tRCS)       {O = tQ; R =tR;}
//    CoordSys(CoordSys *tRCS, const Matrix3 &tO, const Vector3 &tR)  : CoordSys(tRCS)    {O = tO; R =tR;}

    //--- Getters

    CoordSys *Get_Ref_CoordSys()    {return Ref_CoordSys;}

    //--- Positions

    Vector3 Position_Local()            {return R;}
    Vector3 Position_Global()
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Position_Global(R);
        else                return R;
    }
    Vector3 Position_Global(const Vector3 &RRel)
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Position_Global(Orient_Mat()*RRel + R);
        else                return RRel;
    }

    //--- Orientations

    Quat Orientation_Global()
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Orientation_Global(O);
        else                return O;
    }

    Quat Orientation_Global(const Quat &ORel)
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Orientation_Global(O*ORel*O.inverse());
        else                return O*ORel*O.inverse();
    }

    Matrix3 OrientMat(const Quat &Q)
    {
        Matrix3 M;
        M.row(0) = Q*UnitX;
        M.row(1) = Q*UnitY;
        M.row(2) = Q*UnitZ;
        return M;
    }

    void Set_Local_Orientation(Vector3 &nX, Vector3 &nY, Vector3 &nZ)
    {
        Matrix3 M;
        M.col(0) = nX;
        M.col(1) = nY;
        M.col(2) = nZ;
        O = Quat(M);
    }

    void Rotate_Centroid_about_z(Real &th){
        // Small helper function for slightly rotating
        Vector3 x = X_Axis();
        Vector3 y = Y_Axis();
        Vector3 xn =  cos(th)*x + sin(th)*y;
        Vector3 yn = -sin(th)*x + cos(th)*y;
        Vector3 zn = xn.cross(yn);
        Set_Local_Orientation(xn,yn,zn);
    }

    Vector3 X_Axis()            {return O*UnitX;}
    Vector3 Y_Axis()            {return O*UnitY;}
    Vector3 Z_Axis()            {return O*UnitZ;}
    Matrix3 OrientMat_Local()   {return OrientMat(O);}

    Vector3 X_Axis_Global()     {return Orientation_Global()*UnitX;}
    Vector3 Y_Axis_Global()     {return Orientation_Global()*UnitY;}
    Vector3 Z_Axis_Global()     {return Orientation_Global()*UnitZ;}
    Matrix3 OrientMat_Global()  {return OrientMat(Orientation_Global());}

    //--- Velocities

    Vector3 Velocity_Local()            {return U;}
    Vector3 Velocity_Global()
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Velocity_Global(R,U);
        else                return U;
    }

    Vector3 Velocity_Global(const Vector3 &RRel, const Vector3 &URel)
    {
        Vector3 U_W_Rel = URel + U_Rot_Rel(RRel) + U;
        if (Ref_CoordSys)   return Ref_CoordSys->Velocity_Global(R,Orient_Mat()*U_W_Rel);
        else                return U_W_Rel;
    }

    Vector3 U_Rot_Rel(const Vector3 &RRel)
    {
        Matrix3 RM;
        RM << 0, -RRel(2), RRel(1), RRel(2), 0, -RRel(0), -RRel(1), RRel(0), 0;
        return -RM*W;
    }

    //--- Angular velocities

    Matrix3 Rel_Pos_Mat()       {Matrix3 RM; RM << 0, -R(2), R(1), R(2), 0, -R(0), -R(1), R(0), 0; return RM;}
    Vector3 Angular_Velocity_Local()    {return W;}
    Vector3 Angular_Velocity_Global()
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Angular_Velocity_Global(W);
        else                return U;
    }
    Vector3 Angular_Velocity_Global(const Vector3 &WRel)
    {
        if (Ref_CoordSys)   return Ref_CoordSys->Angular_Velocity_Global(Orient_Mat()*WRel + W);
        else                return WRel;
    }

    //--- Helper Fn's

    Matrix3 Orient_Mat()                {return O.toRotationMatrix();}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}

#endif // COORDSYS_H
