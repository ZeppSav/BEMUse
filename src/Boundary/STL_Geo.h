/****************************************************************************
    BEMUse Boundary
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

    -> This contains the class which imports the geometry from an stl file.

*****************************************************************************/

#ifndef STL_GEO_H
#define STL_GEO_H

#include "Boundary_Base.h"

namespace BEMUse
{

class STL_Geometry : public Boundary
{
protected:

    std::vector<std::vector<Vector3>> STL_Els_Raw;
    std::vector<std::vector<SP_Node>> STL_Els;

public:

    //--- Constructor
    STL_Geometry()  {}

    //--- Import/Export
    bool Read_Input_File(std::string &FilePath);

    //--- Geometry generation
    void Generate_Nodes();
//    void Generate_Aux_Nodes()       {}
    void Generate_Elements();
//    void Generate_Aux_Elements()    {}

};

}

#endif // STL_GEO_H
