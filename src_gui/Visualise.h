/****************************************************************************
    BEMUse Visualise
    Copyright (C) 2022 Joseph Saverin j.saverin@tu-berlin.de

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

    Information on file:

    -> Visualisations functions for boundary elements

*****************************************************************************/

#ifndef VISUALISE_H
#define VISUALISE_H

#include "src/Boundary/Boundary_Base.h"
#include "src/Geometry/Panel.h"

#ifdef _WIN32
#include <windows.h>
#endif

#include <GL/glu.h>

namespace BEMUse
{

    static GLfloat LineWidth = 0.5;
    static GLfloat LineOpacity = 1.0;
    static GLfloat PanelOpacity = 0.75;
    static GLfloat WaveLineOpacity = 0.25;
    static GLfloat WavePanelOpacity = 0.25;

    static void Vis_Bi_Outline(StateVector &S, Real Op)
    {
        // Vis 2D element
        glColor4d(0,0,0,Op);         // Surface outlines are black
        glLineWidth(LineWidth);
        glBegin(GL_LINES);
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glEnd();
    }

    static void Vis_Tri_Outline(StateVector &S, Real Op)
    {
        // Vis 3D tri element
        glColor4d(0,0,0,Op);         // Surface outlines are black
        glLineWidth(LineWidth);
        glBegin(GL_LINES);
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glVertex3d(S[2](0),S[2](1),S[2](2));
        glVertex3d(S[2](0),S[2](1),S[2](2));
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glEnd();
    }

    static void Vis_Quad_Outline(StateVector &S, Real Op)
    {
        // Vis 3D quad element
        glColor4d(0,0,0,Op);         // Surface outlines are black
        glLineWidth(LineWidth);
        glBegin(GL_LINES);
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glVertex3d(S[2](0),S[2](1),S[2](2));
        glVertex3d(S[2](0),S[2](1),S[2](2));
        glVertex3d(S[3](0),S[3](1),S[3](2));
        glVertex3d(S[3](0),S[3](1),S[3](2));
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glEnd();
    }

    static void Vis_Outline(StateVector &S, Real Op=LineOpacity)
    {
        if (S.size()==2)    Vis_Bi_Outline(S,Op);
        if (S.size()==3)    Vis_Tri_Outline(S,Op);
        if (S.size()==4)    Vis_Quad_Outline(S,Op);
    }

    static void Vis_Outline(SP_Geo E, Real Op=LineOpacity)
    {
        StateVector S;
        E->Get_Nodal_Positions(S);
        Vis_Outline(S,Op);
    }

    static void Visualise_Boundary_Elements(Boundary *B)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Els;
        B->Get_Elements(Els);
//        B->Get_Symm_Elements(Els);
        for (SP_Geo E : Els) Vis_Outline(E);
    }

    static void Visualise_Interior_FreeSurface_Elements(Boundary *B)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Els;
        B->Get_Aux_Elements(Els);
        for (SP_Geo E : Els) Vis_Outline(E);
    }

    static void Visualise_Exterior_FreeSurface_Elements(Boundary *B)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Els;
        B->Get_FreeSurface_Elements(Els);
        for (SP_Geo E : Els) Vis_Outline(E);
    }

    static void Vis_CS(SP_Node N)
    {
        // This simply visualises the elements
        Vector3 P1 = N->Position_Global();
//        Vector3 P2 = N->X_Axis_Global();
//        Vector3 P2 = N->Y_Axis_Global();
//        Vector3 P2 = N->Z_Axis_Global();
        Vector3 P2 = N->Z_Axis_Global();

        glColor4d(0,0,0,LineOpacity);         // Surface outlines are black
        glLineWidth(LineWidth);
        Real L = 0.1;

        // Plot vector
        glBegin(GL_LINES);
        glVertex3d(P1(0),P1(1),P1(2));
        glVertex3d(P1(0)+L*P2(0),P1(1)+L*P2(1),P1(2)+L*P2(2));
        glEnd();
    }

    static void Visualise_Centroid_CS(Boundary *B)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Gs;
        B->Get_Elements(Gs);
//        B->Get_Symm_Elements(Gs);
        for (SP_Geo G : Gs) G->Set_Centroid();
        for (SP_Geo G : Gs) Vis_CS(G->Centroid);
    }

    static void Visualise_Node_CS(Boundary *B)
    {
        // This simply visualises the elements
        std::vector<SP_Node> Ns;
        B->Get_Nodes(Ns);
        for (SP_Node N : Ns) Vis_CS(N);
    }

    static void Vis_Bi_Panel(StateVector &S, Vector3 &C, Real Op)
    {
        // Vis 2D element
        glColor4d(C(0),C(1),C(2),Op);         // Surface outlines are black
        glLineWidth(LineWidth);
        glBegin(GL_LINES);
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glEnd();
    }

    static void Vis_Tri_Panel(StateVector &S, Vector3 &C, Real Op)
    {
        // Vis 3D tri element
        glColor4d(C(0),C(1),C(2),Op);         // Surface outlines are black
        glBegin(GL_TRIANGLES);
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glVertex3d(S[2](0),S[2](1),S[2](2));
        glEnd();
    }

    static void Vis_Quad_Panel(StateVector &S, Vector3 &C, Real Op)
    {
        // Vis 3D quad element
        glColor4d(C(0),C(1),C(2),Op);
        glBegin(GL_QUADS);
        glVertex3d(S[0](0),S[0](1),S[0](2));
        glVertex3d(S[1](0),S[1](1),S[1](2));
        glVertex3d(S[2](0),S[2](1),S[2](2));
        glVertex3d(S[3](0),S[3](1),S[3](2));
        glEnd();
    }

    static void Vis_Panel(StateVector S, Vector3 C, Real Op=PanelOpacity)
    {
        if (S.size()==2)    Vis_Bi_Panel(S,C,Op);
        if (S.size()==3)    Vis_Tri_Panel(S,C,Op);
        if (S.size()==4)    Vis_Quad_Panel(S,C,Op);
    }

    static void Vis_Panel(SP_Geo E, Vector3 C)
    {
        StateVector S;
        E->Get_Nodal_Positions(S);
        Vis_Panel(S,C);
    }

    static void Visualise_Solution(Boundary *B)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Els;
        B->Get_Elements(Els);
        for (SP_Geo E : Els) Vis_Panel(E,E->Centroid->VWeight);
    }

    static Vector3 Set_Colour(Real &SMin, Real &SMax, Real &f)
    {
        // Specifies a colour spectrum
        Real val = (f-SMin)/(SMax-SMin);

        // Red-Green scheme
        Vector3 Blue(0,0,1), Red(1,0,0),  Green(0,1,0);
        if (val<0.5)    return Red + 2*val*(Green-Red);
        else            return Green + 2*(val-0.5)*(Blue-Green);
    }

    static void Visualise_Radiation_Solution(Boundary *B, int Freq, int DOF, bool R = true)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Els;
        B->Get_Elements(Els);

        // Prepares surface array
        Vector A;
        if (R)  A = B->RadSolArray[Freq].col(DOF).real();
        else    A = B->RadSolArray[Freq].col(DOF).imag();
        Real AMax = A.maxCoeff(), AMin = A.minCoeff();
        for (int i=0; i<Els.size(); i++) Vis_Panel(Els[i],Set_Colour(AMin,AMax,A[i]));
    }

    static void Visualise_Diffraction_Solution(Boundary *B, int Freq, int Beta, bool R = true)
    {
        // This simply visualises the elements
        std::vector<SP_Geo> Els;
        B->Get_Elements(Els);

        // Prepares surface array
        Vector A;
        if (R)  A = B->DiffSolArray[Freq].col(Beta).real();
        else    A = B->DiffSolArray[Freq].col(Beta).imag();
        Real AMax = A.maxCoeff(), AMin = A.minCoeff();
        for (int i=0; i<Els.size(); i++) Vis_Panel(Els[i],Set_Colour(AMin,AMax,A[i]));
    }

    static void Visualise_Free_Surface(Boundary *B, int Freq, int DOF, Real Time)
    {
        // This visualises the free surface
        // Collect node positions
        std::vector<SP_Node> Nodes;
        B->Get_FreeSurface_Nodes(Nodes);
        if (Nodes.empty())    return;

        Real TFac = (exp(CReal(0,Time))).real();

        StateVector FSNodePos;
        for (int i=0; i<Nodes.size(); i++){
            Vector3 FSPos = Nodes[i]->Position_Global();
            CReal Psi = B->FS_Rad_Array[Freq](Nodes[i]->ID,DOF);
//            CReal Psi = B->FS_Scat_Array[Freq](Nodes[i]->ID,DOF);
            FSNodePos.push_back(FSPos + Vector3(0,0,TFac*Psi.real()));
        }

        // Vis Panels
        std::vector<SP_Geo> Els;
        B->Get_FreeSurface_Elements(Els);
        if (Els.empty())    return;

        for (int i=0; i<Els.size(); i++){
            StateVector Vertices;
            for (int j=0; j<Els[i]->Get_N(); j++) Vertices.push_back(FSNodePos[Els[i]->Get_Node(j)->ID]);
            Vis_Panel(Vertices,Vector3(0,0,1),WavePanelOpacity);
            Vis_Outline(Vertices,WaveLineOpacity);
        }
    }

    static void Visualise_ScatteredFree_Surface(Boundary *B, int Freq, int Beta, Real Time)
    {
        // This visualises the free surface
        // Collect node positions
        std::vector<SP_Node> Nodes;
        B->Get_FreeSurface_Nodes(Nodes);
        if (Nodes.empty())    return;

        Real TFac = (exp(CReal(0,Time))).real();

        StateVector FSNodePos;
        for (int i=0; i<Nodes.size(); i++){
            Vector3 FSPos = Nodes[i]->Position_Global();
//            CReal Psi = B->FS_Rad_Array[Freq](Nodes[i]->ID,DOF);
            CReal Psi = B->FS_Scat_Array[Freq](Nodes[i]->ID,Beta);
            FSNodePos.push_back(FSPos + Vector3(0,0,TFac*Psi.real()));
        }

        // Vis Panels
        std::vector<SP_Geo> Els;
        B->Get_FreeSurface_Elements(Els);
        if (Els.empty())    return;

        for (int i=0; i<Els.size(); i++){
            StateVector Vertices;
            for (int j=0; j<Els[i]->Get_N(); j++) Vertices.push_back(FSNodePos[Els[i]->Get_Node(j)->ID]);
            Vis_Panel(Vertices,Vector3(0,0,1),WavePanelOpacity);
            Vis_Outline(Vertices,WaveLineOpacity);
        }
    }

    static void Heatmap(CMatrix &M)
    {
        // This function plots a matrix to demonstrate where the largest entries are
        // Simply a test function to see what the distribution in the matrices look like

        int nR = M.rows(), nC = M.cols();
        Real max = M.cwiseAbs().maxCoeff(), min = M.cwiseAbs().minCoeff();

        // Plot sections
        for (int i=0; i<nR; i++){
            Real x1 = -2.0+4.0*i/nR;
            Real x2 = -2.0+4.0*(i+1)/nR;
            for (int j=0; j<nC; j++){
                Real y1 = -2.0+4.0*j/nC;
                Real y2 = -2.0+4.0*(j+1)/nC;
                StateVector S;
                S.push_back(Vector3(x1,y1,0.0));
                S.push_back(Vector3(x1,y2,0.0));
                S.push_back(Vector3(x2,y2,0.0));
                S.push_back(Vector3(x2,y1,0.0));
                Real smag = std::abs(M(i,j));
                Vis_Panel(S,Set_Colour(min,max,smag));
            }
        }


    }
}

#endif // VISUALISE_H
