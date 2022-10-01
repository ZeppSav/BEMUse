/****************************************************************************
    BEMUse Library
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

    -> This is entry point for the BEMUse solver

*****************************************************************************/

#ifdef BEMUse_GUI

//----------- BEMUse: Compilation of GUI/Command-line version with qmake ------------

#include "src_gui/BEMUser_Interface.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    BEMUser_Interface w;
    w.show();
    return a.exec();
}

#else

//----------- BEMUse: Compilation of Command-line version with CMake ------------

#include "src/BEMUse_Inquiry.h"
#include "src/BEMUser_Console.h"

int main(int argc, char *argv[])
{
    // Step 1: Is this simply an inquiry for information? If so, avoid following steps.
    if (BEMUse_Enquiry(argc,argv)) return 0;

    // Step 2: This is not an enquiry... proceed with analysis.
    BEMUse_Console C;
    C.Specify_Compute_Params(argc,argv);
    C.Specify_Geometry(argc,argv);
    C.Specify_Solver(argc,argv);
    C.Execute();
}

#endif
