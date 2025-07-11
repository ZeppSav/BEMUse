//-----------------------------------------------------------------------------
//-------------------------Geometry Element Functions--------------------------
//-----------------------------------------------------------------------------

#include "Geo_Elements.h"
#include "Quadrule.cpp"

namespace BEMUse
{

//--- Base class

Real Geometry_Element::Get_Lmin()
{
    Real L = 1e6;
    int N = Get_N();
    StateVector P;
    for (int i=0; i<N; i++) P.push_back(Nodes[i]->Position_Global());
    for (int i=0; i<N; i++){
        for (int j=i; j<N; j++){
            Vector3 D = P[i]-P[j];
            if (D.norm()<L) L =  D.norm();
        }
    }
    return L;
}

Real Geometry_Element::Get_Lmax()
{
    Real L = -1e6;
    int N = Get_N();
    StateVector P;
    for (int i=0; i<N; i++) P.push_back(Nodes[i]->Position_Global());
    for (int i=0; i<N; i++){
        for (int j=i; j<N; j++){
            Vector3 D = P[i]-P[j];
            if (D.norm()>L) L = D.norm();
        }
    }
    return L;
}

//--- Triangular elements

void Tri_Element::Set_Centroid()
{
    // Creates the centroid of the panel
    if (Centroid!=nullptr) return;

    // Takes the mean position of the three nodes
    Vector3 P0 = Nodes[0]->Position_Local();
    Vector3 P1 = Nodes[1]->Position_Local();
    Vector3 P2 = Nodes[2]->Position_Local();
    Vector3 PC = (P0+P1+P2)/3.0;

    Vector3 X1 = P0-0.5*(P1+P2);
    Vector3 X2 = P2-P1;
    Vector3 X3 = X1.cross(X2);

    Area = 0.5*X3.norm();
//    X1.normalize();
//    X2.normalize();
    X3.normalize();

    CoordSys *CSRef = Nodes[0]->Get_Ref_CoordSys();
    Centroid = std::make_shared<Node>(CSRef,Quat::FromTwoVectors(UnitZ,X3),PC);
}

void Tri_Element::Set_Quad_Nodes()
{
    // This function specifies the position and strengths of the quadrature nodes.
    if (!QNodes.empty()) return;

    // Gauss Legendre (Dunavant)
    int N_Points = dunavant_order_num(N);
    double xy[2*N_Points], W[N_Points],  PX[N_Points],PY[N_Points];
    dunavant_rule(N, N_Points, xy, W);
    for(int i=0; i<N_Points; i++){
        PX[i] = xy[i*2];
        PY[i] = xy[i*2+1];
        W[i] *= 0.5;
    }

    // Create nodes with the given interpolated coordinate systems.
    for (int i=0; i<N_Points; i++)
    {
        //Nodal Shape Functions for Triangle

        Real xi   = PX[i];
        Real eta  = PY[i];

        Vector3 P(0,0,0);
        P += Nodes[0]->Position_Local()*(1.-xi-eta);
        P += Nodes[1]->Position_Local()*(eta);
        P += Nodes[2]->Position_Local()*(xi);

        // Surface normals
        Vector3 R_xi = Vector3::Zero() , R_eta = Vector3::Zero();
        R_xi += Nodes[0]->Position_Local()*(-1);
        R_xi += Nodes[2]->Position_Local()*(1);
        R_eta += Nodes[0]->Position_Local()*(-1);
        R_eta += Nodes[1]->Position_Local()*(1);

        Vector3 R_n = R_xi.cross(R_eta);
        Real Jac = R_n.norm();              // Jacobian of quadrature node. Must be multiplied with weight.

//        R_xi.normalize();
//        R_eta.normalize();
        R_n.normalize();

        SP_Node QNode = std::make_shared<Node>(Centroid->Get_Ref_CoordSys(),Quat::FromTwoVectors(UnitZ,R_n),P);
        QNode->Weight = Jac*W[i];
        QNodes.push_back(QNode);

        // Store quadrature positions
        QNode->QuadPos(0) = PX[i];
        QNode->QuadPos(1) = PY[i];
    }
}

//--- Quadrilateral elements

void Quad_Element::Set_Centroid()
{
    // Creates the centroid of the panel
    if (Centroid!=nullptr) return;

    // Takes the mean position of the four nodes
    Vector3 P0 = Nodes[0]->Position_Local();
    Vector3 P1 = Nodes[1]->Position_Local();
    Vector3 P2 = Nodes[2]->Position_Local();
    Vector3 P3 = Nodes[3]->Position_Local();
    Vector3 PC = 0.25*(P0+P1+P2+P3);

    Vector3 X1 = 0.5*(P3+P2)-0.5*(P0+P1);       // X1.normalize();
    Vector3 X2 = 0.5*(P2+P1)-0.5*(P3+P0);       // X2.normalize();
    Vector3 X3 = X1.cross(X2);                  X3.normalize();

    Vector3 C1 = P2-P0;
    Vector3 C2 = P3-P1;
    Area = 0.5*(C1.cross(C2)).norm();

    CoordSys *CSRef = Nodes[0]->Get_Ref_CoordSys();
    Centroid = std::make_shared<Node>(CSRef,Quat::FromTwoVectors(UnitZ,X3),PC);
}

void Quad_Element::Set_Quad_Nodes()
{
    // This function specifies the position and strengths of the quadrature nodes.
    if (!QNodes.empty()) return;

    // Depends on qadrature discretisation.
    Cart_ID2 QS = Set_Quad_Shape();

    double WX[QS(0)],PX[QS(0)],WY[QS(1)],PY[QS(1)];

    // Gauss Legendre
    cdgqf(QS(0),1,0,0,PX,WX);
    cdgqf(QS(1),1,0,0,PY,WY);

    // Chebyshev quad
//    chebyshev_set(QS(0),PX,WX);
//    chebyshev_set(QS(1),PY,WY);

    // Gauss Lobatto
//    lobatto_compute(QS(0),PX,WX);
//    lobatto_compute(QS(1),PY,WY);

    // Clenshaw curtis
//    clenshaw_curtis_compute(QS(0),PX,WX);
//    clenshaw_curtis_compute(QS(1),PY,WY);

//    for (int i=0; i<QS(0); i++) qDebug() << PX[i] << WX[i];

    // Create nodes with the given interpolated coordinate systems.
    for (int i=0; i<QS(0); i++){
        for (int j=0; j<QS(1); j++){

            Vector3 P(0,0,0);
            P += Nodes[0]->Position_Local()*0.25*(1-PX[i])*(1-PY[j]);
            P += Nodes[1]->Position_Local()*0.25*(1-PX[i])*(1+PY[j]);
            P += Nodes[2]->Position_Local()*0.25*(1+PX[i])*(1+PY[j]);
            P += Nodes[3]->Position_Local()*0.25*(1+PX[i])*(1-PY[j]);

            // Surface normals
            Vector3 R_xi = Vector3::Zero();
            R_xi += Nodes[0]->Position_Local()*(-0.25*(1-PY[j]));
            R_xi += Nodes[1]->Position_Local()*(-0.25*(1+PY[j]));
            R_xi += Nodes[2]->Position_Local()*(0.25*(1+PY[j]));
            R_xi += Nodes[3]->Position_Local()*(0.25*(1-PY[j]));

            Vector3 R_eta = Vector3::Zero();
            R_eta += Nodes[0]->Position_Local()*(-0.25*(1-PX[i]));
            R_eta += Nodes[1]->Position_Local()*(0.25*(1-PX[i]));
            R_eta += Nodes[2]->Position_Local()*(0.25*(1+PX[i]));
            R_eta += Nodes[3]->Position_Local()*(-0.25*(1+PX[i]));

            Vector3 R_n = R_xi.cross(R_eta);
            Real Jac = R_n.norm();              // Jacobian of quadrature node. Must be multiplied with weight.

//            R_xi.normalize();
//            R_eta.normalize();
            R_n.normalize();

            SP_Node QNode = std::make_shared<Node>(Centroid->Get_Ref_CoordSys(),Quat::FromTwoVectors(UnitZ,R_n),P);
            QNode->Weight = Jac*WX[i]*WY[j];
            QNodes.push_back(QNode);

            // Store quadrature positions
            QNode->QuadPos(0) = PX[i];
            QNode->QuadPos(1) = PY[j];
        }
    }
}

Cart_ID2 Quad_Element::Set_Quad_Shape()
{
    // This function checks the geometry and determines the optimal quadrature shape/order.
    // For now hard-coded to return simply the desired Quad order
//    int N = 1;
    int NQX = N;
    int NQY = N;

    return Cart_ID2(NQX,NQY);
}

}
