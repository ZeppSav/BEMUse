/****************************************************************************
    Node Class
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

    -> This is the implementation of the node type, which inherits from the coordinate system

*****************************************************************************/


#ifndef NODE_H
#define NODE_H

#include "CoordSys.h"

namespace BEMUse
{

class Node : public CoordSys
{
public:

    //--- Constructors

    Node()                                                      : CoordSys()            {}
    Node(CoordSys *tRCS)                                        : CoordSys(tRCS)        {}
    Node(CoordSys *tRCS, const Vector3 &tR)                     : CoordSys(tRCS,tR)     {}
    Node(CoordSys *tRCS, const Quat &tQ, const Vector3 &tR)     : CoordSys(tRCS,tQ,tR)  {}
//    Node(CoordSys *tRCS, const Matrix3 &tO, const Vector3 &tR)  : CoordSys(tRCS,tO,tR)  {}

    Real    Weight = 0;
    CReal   CWeight = CReal(0,0);
    Vector3 VWeight = Vector3::Ones();
    unsigned int    ID = 0;
};

typedef std::shared_ptr<Node> SP_Node;
}

#endif // NODE_H
