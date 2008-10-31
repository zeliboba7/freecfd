/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#ifndef STATE_CACHE_H
#define STATE_CACHE_H

#include "vec3d.h"

class Cell_State {
	public:
		double rho,p,a,H,k,omega,mu,k_center,omega_center;
		Vec3D v,v_center,vN;
};

class Face_State {
	public:
		double rho,p,k,omega,mu;
		Vec3D v;
		Vec3D gradU,gradV,gradW,gradK,gradOmega;
		Vec3D normal,tangent1,tangent2,left2right;
		double area;
		int bc;
};

#endif
