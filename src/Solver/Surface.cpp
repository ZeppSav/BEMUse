//------------------------------------------------------------
//------------------- Surface Functions ----------------------
//------------------------------------------------------------

#include "Surface.h"

namespace BEMUse
{

//--- Surface class

Surface::Surface(std::vector<SP_Geo> &SurfPans, std::vector<SP_Node> &SurfNodes, DistType Type)
{
    // Constructor to define the surface
    PanDist = Type;

    // Generate panels
    for (SP_Geo G : SurfPans){
        if (G->Get_Type()==TRI)  Panels.push_back(std::make_shared<FlatSourceTriPanel>(G));
        if (G->Get_Type()==QUAD) Panels.push_back(std::make_shared<FlatSourceQuadPanel>(G));
        G->ID = NP++;
    }

    // Store nodes
    StdAppend(Nodes,SurfNodes);
    NN = size(Nodes);

    // Specify panel Connectivity for panel gradient approach
    PG_Conn.resize(NN);
    PG_Angles.resize(NN);
    for (int i=0; i<NP; i++){
        if (Panels[i]->Get_Type()==TRI)  Specify_Connectivity_Tri_Panel(Panels[i]);
        if (Panels[i]->Get_Type()==QUAD) Specify_Connectivity_Quad_Panel(Panels[i]);
    }
    std::cout << "Surface::Specify_Connectivity_Panel_Gradient(): Node-Panel connectivity has been specified for panel gradient approach" << std::endl;

    // Specify panel connectivity for B-spline approach
    Specify_Connectivity_BSpline();

    // Initialize pEta array
    purtEta = Matrix(NN,1);
}

// Surface distribution: Panel gradient approach

void Surface::Specify_Connectivity_Tri_Panel(SP_Panel G)
{
    SP_Node ND1 = G->Get_Geo()->Get_Node(0);
    SP_Node ND2 = G->Get_Geo()->Get_Node(1);
    SP_Node ND3 = G->Get_Geo()->Get_Node(2);

    Vector3 Pos1 = ND1->Position_Global();
    Vector3 Pos2 = ND2->Position_Global();
    Vector3 Pos3 = ND3->Position_Global();

    // Add this panel to this node's list of Connected panels.
    PG_Conn[ND1->ID].push_back(G->Get_Geo()->ID);
    PG_Conn[ND2->ID].push_back(G->Get_Geo()->ID);
    PG_Conn[ND3->ID].push_back(G->Get_Geo()->ID);

    // Calculate connection angles.
    Vector3 nd1a = Pos2-Pos1, nd1b = Pos3-Pos1;
    Vector3 nd2a = Pos1-Pos2, nd2b = Pos3-Pos2;
    Vector3 nd3a = Pos1-Pos3, nd3b = Pos2-Pos3;
    PG_Angles[ND1->ID].push_back(acos(nd1a.dot(nd1b)/(nd1a.norm()*nd1b.norm())));
    PG_Angles[ND2->ID].push_back(acos(nd2a.dot(nd2b)/(nd2a.norm()*nd2b.norm())));
    PG_Angles[ND3->ID].push_back(acos(nd3a.dot(nd3b)/(nd3a.norm()*nd3b.norm())));
}

void Surface::Specify_Connectivity_Quad_Panel(SP_Panel G)
{
    SP_Node ND1 = G->Get_Geo()->Get_Node(0);
    SP_Node ND2 = G->Get_Geo()->Get_Node(1);
    SP_Node ND3 = G->Get_Geo()->Get_Node(2);
    SP_Node ND4 = G->Get_Geo()->Get_Node(3);

    Vector3 Pos1 = ND1->Position_Global();
    Vector3 Pos2 = ND2->Position_Global();
    Vector3 Pos3 = ND3->Position_Global();
    Vector3 Pos4 = ND4->Position_Global();

    // Add this panel to this node's list of Connected panels.
    PG_Conn[ND1->ID].push_back(G->Get_Geo()->ID);
    PG_Conn[ND2->ID].push_back(G->Get_Geo()->ID);
    PG_Conn[ND3->ID].push_back(G->Get_Geo()->ID);
    PG_Conn[ND4->ID].push_back(G->Get_Geo()->ID);

    // Calculate connection angles.
    Vector3 nd1a = Pos2-Pos1, nd1b = Pos4-Pos1;
    Vector3 nd2a = Pos1-Pos2, nd2b = Pos3-Pos2;
    Vector3 nd3a = Pos2-Pos3, nd3b = Pos4-Pos3;
    Vector3 nd4a = Pos1-Pos4, nd4b = Pos3-Pos4;
    PG_Angles[ND1->ID].push_back(acos(nd1a.dot(nd1b)/(nd1a.norm()*nd1b.norm())));
    PG_Angles[ND2->ID].push_back(acos(nd2a.dot(nd2b)/(nd2a.norm()*nd2b.norm())));
    PG_Angles[ND3->ID].push_back(acos(nd3a.dot(nd3b)/(nd3a.norm()*nd3b.norm())));
    PG_Angles[ND4->ID].push_back(acos(nd4a.dot(nd4b)/(nd4a.norm()*nd4b.norm())));
}

void Surface::Calculate_PG_Coeffs_Tri_Panel(SP_Panel G)
{
    SP_Node ND1 = G->Get_Geo()->Get_Node(0);
    SP_Node ND2 = G->Get_Geo()->Get_Node(1);
    SP_Node ND3 = G->Get_Geo()->Get_Node(2);

    Vector3 CLoc = G->Get_Geo()->Centroid->Position_Global();
    Matrix3 OLoc = G->Get_Geo()->Centroid->OrientMat_Global();
    Vector3 P_Loc1 = OLoc * (ND1->Position_Global() - CLoc);
    Vector3 P_Loc2 = OLoc * (ND2->Position_Global() - CLoc);
    Vector3 P_Loc3 = OLoc * (ND3->Position_Global() - CLoc);

    Matrix3 Surface_Dist = Matrix3::Ones();
    Surface_Dist(0,1) = P_Loc1(0);    Surface_Dist(0,2) = P_Loc1(1);
    Surface_Dist(1,1) = P_Loc2(0);    Surface_Dist(1,2) = P_Loc2(1);
    Surface_Dist(2,1) = P_Loc3(0);    Surface_Dist(2,2) = P_Loc3(1);
    PG_Mats[G->Get_Geo()->ID] = Surface_Dist.inverse();
}

void Surface::Calculate_PG_Coeffs_Quad_Panel(SP_Panel G)
{
    SP_Node ND1, ND2, ND3, ND4;
    ND1 = G->Get_Geo()->Get_Node(0);
    ND2 = G->Get_Geo()->Get_Node(1);
    ND3 = G->Get_Geo()->Get_Node(2);
    ND4 = G->Get_Geo()->Get_Node(3);

    Vector3 CLoc = G->Get_Geo()->Centroid->Position_Global();
    Matrix3 OLoc = G->Get_Geo()->Centroid->OrientMat_Global();
    Vector3 P_Loc0 = OLoc * (ND1->Position_Global() - CLoc);
    Vector3 P_Loc1 = OLoc * (ND2->Position_Global() - CLoc);
    Vector3 P_Loc2 = OLoc * (ND3->Position_Global() - CLoc);
    Vector3 P_Loc3 = OLoc * (ND4->Position_Global() - CLoc);

    Matrix Surface_Dist = Matrix::Ones(4,4);
    Surface_Dist(0,1) = P_Loc0(0);    Surface_Dist(0,2) = P_Loc0(1);    Surface_Dist(0,3) = P_Loc0(0)*P_Loc0(1);
    Surface_Dist(1,1) = P_Loc1(0);    Surface_Dist(1,2) = P_Loc1(1);    Surface_Dist(1,3) = P_Loc1(0)*P_Loc1(1);
    Surface_Dist(2,1) = P_Loc2(0);    Surface_Dist(2,2) = P_Loc2(1);    Surface_Dist(2,3) = P_Loc2(0)*P_Loc2(1);
    Surface_Dist(3,1) = P_Loc3(0);    Surface_Dist(3,2) = P_Loc3(1);    Surface_Dist(3,3) = P_Loc3(0)*P_Loc3(1);
    PG_Mats[G->Get_Geo()->ID] = Surface_Dist.inverse();
}

void Surface::Calculate_PG_Coeffs()
{
    // The coefficients are calculated for each panel which allow for the interpolation of the surface distribution
    PG_Mats.clear();
    PG_Mats.resize(NP);
    OpenMPfor
    for (int i=0; i<NP; i++){
        if (Panels[i]->Get_Type()==TRI)   Calculate_PG_Coeffs_Tri_Panel(Panels[i]);
        if (Panels[i]->Get_Type()==QUAD)  Calculate_PG_Coeffs_Quad_Panel(Panels[i]);
    }
}

void Surface::Calculate_PG_Gradients(Vector &Field, std::vector<Vector3> &Gradients)
{
    // Using the gradients which were specified in the Calculate_PG_Coeffs function, the
    // weighted gradients on the nodes are calculated here.

    // Set the value of the gradient on each panel
    std::vector<Vector3> Panel_Gradients(NP);
    OpenMPfor
    for (int i=0; i<NP; i++){

        SP_Geo G = Panels[i]->Get_Geo();

        // Solve for surface distribution
        if (G->Get_Type()==TRI){         // Triangular panels
            Vector3 SurfDist = Vector3( Field(G->Get_Node(0)->ID),
                                        Field(G->Get_Node(1)->ID),
                                        Field(G->Get_Node(2)->ID));
            Vector3 Sol = PG_Mats[i]*SurfDist;
            Panel_Gradients[i] = Sol(1)*G->Centroid->X_Axis_Global() +
                                 Sol(2)*G->Centroid->Y_Axis_Global();
        }
        if (G->Get_Type()==QUAD){        // Quadrilateral panels
            Vector SurfDist = Vector::Zero(4);
            SurfDist << Field(G->Get_Node(0)->ID),
                        Field(G->Get_Node(1)->ID),
                        Field(G->Get_Node(2)->ID),
                        Field(G->Get_Node(3)->ID);
            Vector Sol = PG_Mats[i]*SurfDist;
            Panel_Gradients[i] = Sol(1)*G->Centroid->X_Axis_Global() +
                                 Sol(2)*G->Centroid->Y_Axis_Global();
        }
    }

    // Calculate the weighted value of gradient for each node using connectivities
    OpenMPfor
    for (int i=0; i<NN; i++){
        Vector3 Grad = Vector3::Zero();
        Real THTOT = 0;
        for (int n=0; n<size(PG_Angles[i]); n++) THTOT += PG_Angles[i][n];        // Theoretically, this should sum to 2Pi
        for (int n=0; n<size(PG_Angles[i]); n++) Grad += PG_Angles[i][n]/THTOT*Panel_Gradients[PG_Conn[i][n]];
        Vector3 xg = Nodes[i]->X_Axis_Global();
        Vector3 yg = Nodes[i]->Y_Axis_Global();
        Real x = Grad.dot(xg);
        Real y = Grad.dot(yg);

        // Return gradient at node
        Gradients[i] = x*xg + y*yg;
    }
}

// Surface distribution: B-Spline approach

void Surface::Specify_Connectivity_BSpline()
{
    // For interpolations of variables over surfaces, the connectivity of the nodes must be specified:
    // This is specified here for the case that the B-Spline approach is used to approximate the local distribution.

    //--- Body Nodes

    // Direct connectivity via panels
    BSP_Conn.resize(NN);
    for (int i=0; i<NP; i++){

        int ID1, ID2, ID3, ID4;
        ID1 = Panels[i]->Get_Geo()->Get_Node(0)->ID;
        ID2 = Panels[i]->Get_Geo()->Get_Node(1)->ID;
        ID3 = Panels[i]->Get_Geo()->Get_Node(2)->ID;
        if (Panels[i]->Get_Type()==QUAD) ID4 = Panels[i]->Get_Geo()->Get_Node(3)->ID;

        if (Panels[i]->Get_Type()==QUAD){
            BSP_Conn[ID1].push_back(ID2);  BSP_Conn[ID1].push_back(ID4);
            BSP_Conn[ID2].push_back(ID1);  BSP_Conn[ID2].push_back(ID3);
            BSP_Conn[ID3].push_back(ID2);  BSP_Conn[ID3].push_back(ID4);
            BSP_Conn[ID4].push_back(ID1);  BSP_Conn[ID4].push_back(ID3);
        }
        if (Panels[i]->Get_Type()==TRI){
            BSP_Conn[ID1].push_back(ID2);  BSP_Conn[ID1].push_back(ID3);
            BSP_Conn[ID2].push_back(ID1);  BSP_Conn[ID2].push_back(ID3);
            BSP_Conn[ID3].push_back(ID1);  BSP_Conn[ID3].push_back(ID2);
        }
    }

    // First neighbors
    OpenMPfor
    for (int i=0; i<NN; i++){
        Sort_Ascending(BSP_Conn[i]);
        Remove_Duplicates(BSP_Conn[i]);
    }

    // Prepare second neigbors (Compose array, sort, remove duplicates, remove first neighbors from second neighbors list)
    BSP_Conn2.resize(NN);
    OpenMPfor
    for (int i=0; i<NN; i++){
        for (int j=0; j<int(BSP_Conn[i].size()); j++){
            int idn1 = BSP_Conn[i][j];                     //  Index of neighbor 1
            StdAppend(BSP_Conn2[i],BSP_Conn[idn1]);
        }
        Sort_Ascending(BSP_Conn2[i]);
        Remove_Duplicates(BSP_Conn2[i]);
        Remove_Value(BSP_Conn2[i],i);
        for (int j=0; j<int(BSP_Conn[i].size()); j++)  Remove_Value(BSP_Conn2[i],BSP_Conn[i][j]);
    }

    // // Debug connected nodes
    // int NNode = NNBD-1;
    // Vector3 pos = Body_Nodes[NNode]->Position_Global();
    // std::cout << "Node " << std::endl;
    // std::cout << pos(0) <<" "<< pos(1) <<" "<< pos(2) << std::endl;
    // std::cout << "First neighbors " << std::endl;
    // for (int j=0; j<int(BSP_Conn[NNode].size()); j++){
    //     pos = Body_Nodes[BSP_Conn[NNode][j]]->Position_Global();
    //     std::cout << pos(0) <<" "<< pos(1) <<" "<< pos(2) << std::endl;
    // }
    // std::cout << "Second neighbors " << std::endl;
    // for (int j=0; j<int(BSP_Conn2[NNode].size()); j++){
    //     pos = Body_Nodes[BSP_Conn2[NNode][j]]->Position_Global();
    //     std::cout << pos(0) <<" "<< pos(1) <<" "<< pos(2) << std::endl;
    // }

    std::cout << "Surface::Specify_Connectivity_BSpline(): Node-node connectivity has been specified for B-Spline approach" << std::endl;
}

void Surface::BSpline_Coefficients()
{
    // This function calculates for each node the matrix of B-Spline coefficients.
    // The matrix corresponding to the calculation of the gradient terms is also calculated

    BSP_Mats.clear();
    BSP_Mats.resize(NN);
    BSP_RBF_Mats.clear();
    BSP_RBF_Mats.resize(NN);

    OpenMPfor
    for (int node=0; node<NN; node++){

        std::vector<SP_Node> all_neighbors;
        all_neighbors.push_back(Nodes[node]);
        for (int n=0; n<size(BSP_Conn[node]); n++)     all_neighbors.push_back(Nodes[BSP_Conn[node][n]]);
        for (int n=0; n<size(BSP_Conn2[node]); n++)    all_neighbors.push_back(Nodes[BSP_Conn2[node][n]]);
        int nrNeigh = all_neighbors.size();
        int system_size = nrNeigh + 6;

        Matrix CoeffMat = Matrix::Zero(system_size, system_size);
        BSP_RBF_Mats[node] = Matrix::Zero(system_size, 2);

        Vector3 xi_u = Nodes[node]->X_Axis_Global();
        Vector3 xi_v = Nodes[node]->Y_Axis_Global();
        // Vector3 xi_Normal = N->Z_Axis_Global();
        Vector3 Cloc = Nodes[node]->Position_Global();

        for (int i = 0; i < nrNeigh; i++) {

            Vector3 xj_i = all_neighbors[i]->Position_Global() - Cloc;

            for (int j = 0; j < nrNeigh; j++) {
                Vector3 xj_j = all_neighbors[j]->Position_Global() - Cloc;
                Vector3 diff = xj_i - xj_j;
                Real diffx = diff.dot(xi_u);
                Real diffy = diff.dot(xi_v);
                Real r = sqrt(diffx * diffx + diffy * diffy);
                CoeffMat(i, j) = r * r * r;
            }

            // Terms for B_Spline coefficients
            Real xloc = xj_i.dot(xi_u);
            Real yloc = xj_i.dot(xi_v);
            CoeffMat(i, nrNeigh + 0) = 1.0;               // Constant term
            CoeffMat(i, nrNeigh + 1) = xloc;              // x term
            CoeffMat(i, nrNeigh + 2) = yloc;              // y term
            CoeffMat(i, nrNeigh + 3) = xloc * yloc;       // xy term
            CoeffMat(i, nrNeigh + 4) = xloc * xloc;       // x^2 term
            CoeffMat(i, nrNeigh + 5) = yloc * yloc;       // y^2 term
            CoeffMat.block(nrNeigh, i, 6, 1) = CoeffMat.block(i, nrNeigh, 1, 6).transpose();

            // Terms for RBF grad coefficients (remember, we are taking xj_i - xj_j)
            Real xrbf = -xloc, yrbf = -yloc, r = sqrt(xrbf*xrbf + yrbf*yrbf);
            BSP_RBF_Mats[node](i,0) = 3.0*xrbf*r;
            BSP_RBF_Mats[node](i,1) = 3.0*yrbf*r;
        }
        BSP_Mats[node] = CoeffMat.inverse();
        // BSP_Mats[node] = CoeffMat;
    }
}

void Surface::BSpline_Gradient(Vector &Field, std::vector<Vector3> &Gradients)
{
    // This function calculates for a given node (N), its neighbors (Ne) and it's second neighbors (Ne2) the matrix
    // of B-Spline coefficients. It also calculates the matrix corresponding to the calculation of the gradient

    OpenMPfor
    for (int node=0; node<NN; node++){

        std::vector<SP_Node> all_neighbors;
        all_neighbors.push_back(Nodes[node]);
        for (int n=0; n<size(BSP_Conn[node]); n++) all_neighbors.push_back(Nodes[BSP_Conn[node][n]]);
        for (int n=0; n<size(BSP_Conn2[node]); n++) all_neighbors.push_back(Nodes[BSP_Conn2[node][n]]);
        int nrNeigh = all_neighbors.size();
        int system_size = nrNeigh + 6;

        // Fill RHS Matrix
        Matrix B = Matrix::Zero(system_size, 1);
        for (int i=0; i<nrNeigh; i++) B(i) = Field(all_neighbors[i]->ID);

        // Solve system
        Matrix coeffs = BSP_Mats[node]*B;

        // We now have the coefficients of the system. We now use these to solve for the gradient (in the local system)
        Real GradXloc = 0., GradYloc = 0.;
        for (int i=0; i<nrNeigh; i++){      // RBF terms
            GradXloc += BSP_RBF_Mats[node](i,0)*coeffs(i);
            GradYloc += BSP_RBF_Mats[node](i,1)*coeffs(i);
        }
        GradXloc += coeffs(nrNeigh+1);      // Polynomial term x
        GradYloc += coeffs(nrNeigh+2);      // Polynomial term y

        // Convert to global system
        Vector3 xl = Nodes[node]->X_Axis_Global();
        Vector3 yl = Nodes[node]->Y_Axis_Global();
        Gradients[node] = GradXloc*xl + GradYloc*yl;
    }
}

// Surface integration

Real Surface::Integrate_Quantities_Quadrature(Matrix &Field, Real exp)
{
    // This function loops over surface elements and uses quadrature to calculate the integral over the surface
    // for triangular elements

    Real I = 0;
    for (SP_Panel Pan : Panels){
        SP_Geo Geo = Pan->Get_Geo();
        for (SP_Node N : Geo->QNodes){
            Real xi = N->QuadPos(0);
            Real eta  = N->QuadPos(1);
            Real t1 = 0.;
            if (Geo->Get_Type()==TRI){                   // Triangular panel
                t1 += pow(Field(Geo->Get_Node(0)->ID),exp)*(1-xi-eta);
                t1 += pow(Field(Geo->Get_Node(1)->ID),exp)*eta;
                t1 += pow(Field(Geo->Get_Node(2)->ID),exp)*xi;
            }
            if (Geo->Get_Type()==QUAD){                   // Quadrilateral panel
                t1 += pow(Field(Geo->Get_Node(0)->ID),exp)*0.25*(1-xi)*(1-eta);
                t1 += pow(Field(Geo->Get_Node(1)->ID),exp)*0.25*(1-xi)*(1+eta);
                t1 += pow(Field(Geo->Get_Node(2)->ID),exp)*0.25*(1+xi)*(1+eta);
                t1 += pow(Field(Geo->Get_Node(3)->ID),exp)*0.25*(1+xi)*(1-eta);
            }
            I += N->Weight*t1;
        }
    }
    return I;
    // std::cout << "Integral over surface has value: " << I << std::endl;
}

// Export visualisation

static int const vtiPrecision = 3;     // Low Precision, low memory
static int const vtiWidth = 11;        // Width for parsing low precision number

void Surface::Export_VTP()
{
    //--- Generate a vtp file of the gemetry
    // For each panel, the individual nod positions are exported and

    // Generate array of nodal positions and panel offsets
    std::vector<Vector3> VNodes;
    std::vector<int> Offsets;
    int NCount = 0;
    for (int i=0; i<NP; i++){
        SP_Geo G = Panels[i]->Get_Geo();
        if (G->Get_Type()==TRI){
            VNodes.push_back(G->Get_Node(0)->Position_Global());
            VNodes.push_back(G->Get_Node(1)->Position_Global());
            VNodes.push_back(G->Get_Node(2)->Position_Global());
            NCount += 3;
        }
        if (G->Get_Type()==QUAD){
            VNodes.push_back(G->Get_Node(0)->Position_Global());
            VNodes.push_back(G->Get_Node(1)->Position_Global());
            VNodes.push_back(G->Get_Node(2)->Position_Global());
            VNodes.push_back(G->Get_Node(3)->Position_Global());
            NCount += 4;
        }
        Offsets.push_back(NCount);
    }

    // Generate array of panel surface strengths
    std::vector<Real> Sigma;
    if (PanDist==CONSTANT){
        for (int i=0; i<NP; i++){
            SP_Geo G = Panels[i]->Get_Geo();
            if (G->Get_Type()==TRI) {
                Sigma.push_back(G->Centroid->VWeight(0));
                Sigma.push_back(G->Centroid->VWeight(0));
                Sigma.push_back(G->Centroid->VWeight(0));
            }
            if (G->Get_Type()==QUAD) {
                Sigma.push_back(G->Centroid->VWeight(0));
                Sigma.push_back(G->Centroid->VWeight(0));
                Sigma.push_back(G->Centroid->VWeight(0));
                Sigma.push_back(G->Centroid->VWeight(0));
            }
        }
    }
    if (PanDist==BILINEAR){
        for (int i=0; i<NP; i++){
            SP_Geo G = Panels[i]->Get_Geo();
            if (G->Get_Type()==TRI) {
                Sigma.push_back(Nodes[G->Get_Node(0)->ID]->VWeight(0));
                Sigma.push_back(Nodes[G->Get_Node(1)->ID]->VWeight(0));
                Sigma.push_back(Nodes[G->Get_Node(2)->ID]->VWeight(0));
            }
            if (G->Get_Type()==QUAD) {
                Sigma.push_back(Nodes[G->Get_Node(0)->ID]->VWeight(0));
                Sigma.push_back(Nodes[G->Get_Node(1)->ID]->VWeight(0));
                Sigma.push_back(Nodes[G->Get_Node(2)->ID]->VWeight(0));
                Sigma.push_back(Nodes[G->Get_Node(3)->ID]->VWeight(0));
            }
        }
    }

    // Export

    CreateDirectory(OutputDirectory);   // Create directory if it doesn't yet exist
    // CreateDirectory(OutputPath);        // Create output file path if it doesn't yet exist

    std::string FilePath = OutputPath + SurfaceName + "_" + std::to_string(currentTimeStep) + ".vtp";
    std::cout << FilePath << std::endl;
    std::ofstream vtifile( FilePath.c_str() );
    vtifile.precision(vtiPrecision);
    if(!vtifile.is_open())
    {
        std::cerr << "ERROR: cannot open vtifile." << std::endl;
        return;
    }

    vtifile << "<?xml version='1.0'?>" << "\n";
    vtifile << "<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>" << "\n";
    vtifile << "  <PolyData>"   << "\n";
    // vtifile << "    <Piece NumberOfPoints='" << NN << "' NumberOfPolys='" << NP << "'> " << "\n";
    vtifile << "    <Piece NumberOfPoints='" << size(VNodes) << "' NumberOfPolys='" << NP << "'> " << "\n";
    vtifile << "        <Points> "   << "\n";
    vtifile << "            <DataArray type='Float32' NumberOfComponents='3' format='ascii'>" << " \n";
    for (size_t i=0; i<size(VNodes); i++){
        // Vector3 P = Nodes[i]->Position_Global();
        // if (fabs(purtEta(i))>100.) purtEta(i) = 0.0;
        // vtifile << std::scientific << std::setw(vtiWidth) << P(0) << std::setw(vtiWidth) << P(1) << std::setw(vtiWidth) << P(2)+purtEta(i);
        vtifile << std::scientific << std::setw(vtiWidth) << VNodes[i](0) << std::setw(vtiWidth) << VNodes[i](1) << std::setw(vtiWidth) << VNodes[i](2);
    }
    vtifile << "            </DataArray> "   << "\n";
    vtifile << "        </Points> "   << "\n";

    vtifile << "        <Polys>"   << "\n";
    vtifile << "            <DataArray type='Int32' Name='connectivity' format='ascii'>" << " \n";
    for (int i=0; i<NCount; i++)    vtifile << std::scientific  << std::setw(vtiWidth) << i;
    vtifile << "            </DataArray> "   << "\n";

    vtifile << "            <DataArray type='Int32' Name='offsets' format='ascii'>" << " \n";
    // for (int i=0; i<NP; i++)  vtifile << std::scientific  << std::setw(vtiWidth) << 3*(i+1);
    for (size_t i=0; i<size(Offsets); i++)  vtifile << std::scientific  << std::setw(vtiWidth) << Offsets[i];
    vtifile << "            </DataArray> "   << "\n";


    vtifile << "        </Polys>"   << "\n";

    vtifile << "            <PointData Scalars='node_scalar'>"   << "\n";
    // if (OutputWaveHeight){
    //     vtifile << "                <DataArray type='Float32' Name='wave_height' format='ascii'>"    << "\n";
    //     for (int i=0; i<NN; i++)  vtifile << std::scientific << std::setw(vtiWidth) << purtEta(i);
    //     vtifile << "                </DataArray>    "       << " \n";
    // }
    vtifile << "                <DataArray type='Float32' Name='node_scalar' format='ascii'>"    << "\n";
    for (size_t i=0; i<size(Sigma); i++)  vtifile << std::scientific  << std::setw(vtiWidth) << Sigma[i];
    vtifile << "                </DataArray>    "       << " \n";
    vtifile << "           </PointData> "   << " \n";

    vtifile << "        </Piece>   "       << " \n";
    vtifile << "    </PolyData>  "       << " \n";
    vtifile << "</VTKFile> "           << " \n";
    vtifile.close();

    std::cout << "Boundary solution (constant panel strength) has been exportet in .vti format to: " << FilePath << std::endl;
}

}
