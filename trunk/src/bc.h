/************************************************************************
	
	Copyright 2007-2010 Emre Sozer

	Contact: emresozer@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#ifndef BC_H
#define BC_H

#include <vector>

// Options for Boundary Condition Types
#define SYMMETRY 1
#define WALL 2
#define INLET 3
#define OUTLET 4
// Options for Boundary Condition Type variants (kind)
#define NONE -1
#define NO_REVERSE 1
#define DAMP_REVERSE 2
#define VELOCITY 3
#define MDOT 4
#define SLIP 5
#define STAGNATION 6
#define FORCE_SUPERSONIC 7
// Options for thermal boundary condition
#define FIXED_T 1
#define FIXED_Q 2
#define ADIABATIC 3
// Options for BC specifications
#define BC_RHO 1
#define BC_P 2
#define BC_T 3
#define BC_STATE 4
#define BC_V 5


class BCregion {
public:
	int type;
	int kind;
	int specified;
	int thermalType;
	double area,total_area;
	Vec3D areaVec;

};

#endif
