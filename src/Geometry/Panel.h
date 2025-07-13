/****************************************************************************
    Panel Class
    Copyright (C) 2021-
    Joseph Saverin j.saverin@tu-berlin.de
    Ronja Baumeister r.baumeister@campus.tu-berlin.de

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

    -> This is the parent class which holds all panel function templates

*****************************************************************************/

#ifndef PANEL_H
#define PANEL_H

#include "Geo_Elements.h"

namespace BEMUse
{

enum DistType {CONSTANT, BILINEAR};

class PanelClass
{
protected:

//    void Set_Quadrature_Nodes();
    SP_Geo Geo;

public:

    //--- Constructors
    PanelClass()     {}
    PanelClass(SP_Geo G)     {Geo=G;}

    //--- Influence constant strength
    virtual Real    Inf_SingleLayer(const Vector3 &P)       {}
    virtual Real    Inf_DoubleLayer(const Vector3 &P)       {}
    virtual void    Inf_SingleDoubleLayer(const Vector3 &P, Real &S, Real &D)   {}

    virtual CReal   CInf_SingleLayer(const Vector3 &P)      {}
    virtual CReal   CInf_DoubleLayer(const Vector3 &P)      {}
    virtual void    CInf_SingleDoubleLayer(const Vector3 &P, CReal &S, CReal &D)   {}

    //--- Influence Gaussian quadrature
    virtual Real    Inf_SingleLayerQuad(const Vector3 &P)   {}
    virtual Real    Inf_DoubleLayerQuad(const Vector3 &P)   {}
    virtual void    Inf_SingleDoubleLayerQuad(const Vector3 &P, Real &S, Real &D)   {}

    virtual CReal   CInf_SingleLayerQuad(const Vector3 &P)  {}
    virtual CReal   CInf_DoubleLayerQuad(const Vector3 &P)  {}
    virtual void    CInf_SingleDoubleLayerQuad(const Vector3 &P, CReal &S, CReal &D)   {}

    //--- Influence non-constant strength
    virtual void Inf_SingleLayer_LinDist(const Vector3 &P, Real &S1, Real &S2, Real &S3) {}
    virtual void Inf_DoubleLayer_LinDist(const Vector3 &P, Real &D1, Real &D2, Real &D3) {}
    virtual void Inf_SingleDoubleLayer_LinDist(const Vector3 &P, Real &S1, Real &S2, Real &S3, Real &D1, Real &D2, Real &D3) {}

    //--- Getters
    SP_Node Centroid()  {return Geo->Centroid;}
    SP_Geo  Get_Geo()   {return Geo;}
    GeoType Get_Type()  {return Geo->Get_Type();}

    //--- Analysis parameters
    Real k = 0.0;       // Frequency parameter

    //--- Specification data
    bool isDipole = false;
};

typedef std::shared_ptr<PanelClass> SP_Panel;

//--- Trilateral panels

class FlatTriPanel : public PanelClass
{
protected:

    //--- Influence
    void    Inf_SourceKernel(const Vector3 &P, Real &S, Real &D);

public:

    //--- Constructors
    FlatTriPanel(SP_Geo G) : PanelClass(G)
    {
        Geo->Set_Centroid();
        Geo->Set_Quad_Nodes();
    }
};

//--- Quadrilateral panels

class FlatQuadPanel : public PanelClass
{
protected:

    //--- Influence
    void    Inf_SourceKernel(const Vector3 &P, Real &S, Real &D);

public:

    //--- Constructors
    FlatQuadPanel(SP_Geo G) : PanelClass(G)
    {
        Geo->Set_Centroid();
        Geo->Set_Quad_Nodes();
    }
};

//--- Source panels

class FlatSourceTriPanel : public FlatTriPanel
{

public:

    //--- Constructors
    FlatSourceTriPanel(SP_Geo G) : FlatTriPanel(G) {}

    //--- Influence constant strength panels
    Real Inf_SingleLayer(const Vector3 &P)   override   {Real S,D; Inf_SourceKernel(P,S,D); return S;}
    Real Inf_DoubleLayer(const Vector3 &P)   override   {Real S,D; Inf_SourceKernel(P,S,D); return D;}
    void Inf_SingleDoubleLayer(const Vector3 &P, Real &S, Real &D) override {Inf_SourceKernel(P,S,D);}

    Real Inf_SingleLayerQuad(const Vector3 &P) override {Real S,D; Inf_SingleDoubleLayerQuad(P,S,D); return S;}
    Real Inf_DoubleLayerQuad(const Vector3 &P) override {Real S,D; Inf_SingleDoubleLayerQuad(P,S,D); return D;}
    void Inf_SingleDoubleLayerQuad(const Vector3 &P, Real &S, Real &D) override;

    //--- Influence linear strength panels
    void Inf_SingleLayer_LinDist(const Vector3 &P, Real &S1, Real &S2, Real &S3) override;
    void Inf_DoubleLayer_LinDist(const Vector3 &P, Real &D1, Real &D2, Real &D3) override;
    void Inf_SingleDoubleLayer_LinDist(const Vector3 &P, Real &S1, Real &S2, Real &S3, Real &D1, Real &D2, Real &D3) override;
};

class FlatSourceQuadPanel : public FlatQuadPanel
{

public:

    //--- Constructors
    FlatSourceQuadPanel(SP_Geo G) : FlatQuadPanel(G) {}

    //--- Influence
    Real Inf_SingleLayer(const Vector3 &P)      {Real S,D; Inf_SourceKernel(P,S,D); return S;}
    Real Inf_DoubleLayer(const Vector3 &P)      {Real S,D; Inf_SourceKernel(P,S,D); return D;}
    void Inf_SingleDoubleLayer(const Vector3 &P, Real &S, Real &D) {Inf_SourceKernel(P,S,D);}

    Real Inf_SingleLayerQuad(const Vector3 &P)  {Real S,D; Inf_SingleDoubleLayerQuad(P,S,D); return S;}
    Real Inf_DoubleLayerQuad(const Vector3 &P)  {Real S,D; Inf_SingleDoubleLayerQuad(P,S,D); return D;}
    void Inf_SingleDoubleLayerQuad(const Vector3 &P, Real &S, Real &D);
};

//--- Wave panels

class FlatWaveTriPanel : public FlatTriPanel
{
protected:

    //--- Testing vars
    Real ZP = -1.0;

public:

    //--- Constructors
    FlatWaveTriPanel(SP_Geo G) : FlatTriPanel(G) {}

    //--- Influence
    CReal CInf_SingleLayer(const Vector3 &P)                                {}
    CReal CInf_DoubleLayer(const Vector3 &P)                                {}
    void  CInf_SingleDoubleLayer(const Vector3 &P, CReal &S, CReal &D)      {}

    CReal CInf_SingleLayerQuad(const Vector3 &P)  {CReal S,D; CInf_SingleDoubleLayerQuad(P,S,D); return S;}
    CReal CInf_DoubleLayerQuad(const Vector3 &P)  {CReal S,D; CInf_SingleDoubleLayerQuad(P,S,D); return D;}
    void  CInf_SingleDoubleLayerQuad(const Vector3 &P, CReal &S, CReal &D);

    void  CInf_SingleDoubleLayerBruteForce(const Vector3 &P, CReal &S, CReal &D);
//    void  Inf_SingleDoubleLayerReflectedSource(const Vector3 &P, CReal &S, CReal &D);
};

class FlatWaveQuadPanel : public FlatQuadPanel
{
protected:

    //--- Testing vars
    Real ZP = -1.0;
//    void Wave_Kernel(const Vector3 &P_glob, const Vector3 &P_Src, const Real &K, CReal &G, CReal &dGX, CReal &dGY, CReal &dGZ);
    void Wave_Kernel_Test();

public:
    //--- Constructors
    FlatWaveQuadPanel(SP_Geo G) : FlatQuadPanel(G) {}

    //--- Influence
    CReal CInf_SingleLayer(const Vector3 &P)                                {}
    CReal CInf_DoubleLayer(const Vector3 &P)                                {}
    void  CInf_SingleDoubleLayer(const Vector3 &P, CReal &S, CReal &D)      {}

    CReal CInf_SingleLayerQuad(const Vector3 &P)  {CReal S,D; CInf_SingleDoubleLayerQuad(P,S,D); return S;}
    CReal CInf_DoubleLayerQuad(const Vector3 &P)  {CReal S,D; CInf_SingleDoubleLayerQuad(P,S,D); return D;}
    void  CInf_SingleDoubleLayerQuad(const Vector3 &P, CReal &S, CReal &D);

    void  CInf_SingleDoubleLayerBruteForce(const Vector3 &P, CReal &S, CReal &D);
//    void  Inf_SingleDoubleLayerSource(const Vector3 &P, CReal &S, CReal &D);
//    void  Inf_SingleDoubleLayerReflectedSource(const Vector3 &P, CReal &S, CReal &D);
};

}

#endif // PANEL_H
