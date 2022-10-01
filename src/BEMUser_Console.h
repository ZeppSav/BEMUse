/****************************************************************************
    BEMUser Console functions
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

    Information on class:

    -> This implements the functionality necesary for running BEMuse from
    -> the command line including input spec

*****************************************************************************/

#ifndef BEMUSER_CONSOLE_H
#define BEMUSER_CONSOLE_H

#include "Boundary/Boundary_Base.h"
#include "Solver/Solver_Base.h"

class BEMUse_Console
{
private:

    //--- BEMUse architecture
    BEMUse::Boundary *Boundary = NULL;
    BEMUse::Solver *Solver = NULL;

    //--- Sim vars
    std::vector<Real> Frequencies;
    std::vector<Real> WaveAngles;

public:

    BEMUse_Console()        {}
    void Specify_Geometry(int argc, char *argv[]);          // Geometry specification
    void Specify_Compute_Params(int argc, char *argv[]);    // Compute parameters
    void Specify_Solver(int argc, char *argv[]);            // Solver specification
    void Execute();                                         // Execute solver
};

#endif // BEMUSER_CONSOLE_H
