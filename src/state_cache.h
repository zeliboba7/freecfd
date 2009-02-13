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
		double p,T,T_center,rho,a,H,k,omega,mu,k_center,omega_center,volume;
		Vec3D v,v_center,vN;
		Vec3D gradU,gradV,gradW,gradK,gradOmega;
		vector<double> update;
		Cell_State &operator= (const Cell_State & rhs) {
			p=rhs.p;
			T=rhs.T;
			T_center=rhs.T_center;
			rho=rhs.rho;
			a=rhs.a;
			H=rhs.H;
			k=rhs.k;
			omega=rhs.omega;
			k_center=rhs.k_center;
			omega_center=rhs.omega_center;
			v=rhs.v;
			gradU=rhs.gradU;
			gradV=rhs.gradV;
			gradW=rhs.gradW;
			gradK=rhs.gradK;
			gradOmega=rhs.gradOmega;
			v_center=rhs.v_center;
			vN=rhs.vN;
			volume=rhs.volume;
			copy(rhs.update.begin(),rhs.update.end(),update.begin());
			return *this;
		}
};

class Face_State {
	public:
		unsigned int index;
		double p,T,rho,k,omega,mu;
		Vec3D v;
		Vec3D gradU,gradV,gradW,gradT,gradK,gradOmega;
		Vec3D normal,tangent1,tangent2,left2right;
		double area;
		int bc;
};

class Fluxes {
	public:
		vector<double> convective,diffusive;
		Fluxes &operator= (const Fluxes & rhs) {
			copy(rhs.convective.begin(),rhs.convective.end(),convective.begin());
			copy(rhs.diffusive.begin(),rhs.diffusive.end(),diffusive.begin());
			return *this;
		}
};

#endif
