/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

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

class BCregion {
public:
	int type;
	int kind;
	int thermalType;
	int specified;
	double T,p,k,omega,rho,mdot,qdot;
	Vec3D v;
	double area;
	Vec3D areaVec;
	double mass,energy;
	Vec3D momentum;
};

class BC {
public:
	vector<BCregion> region;
};

#endif
