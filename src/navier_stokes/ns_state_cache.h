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
#ifndef NS_STATE_CACHE_H
#define NS_STATE_CACHE_H

class NS_Cell_State {
	public:
		double p,p_center,T,T_center,rho,a,H,volume;
		Vec3D V,V_center,Vn;
		vector<double> update;
		NS_Cell_State &operator= (const NS_Cell_State & rhs) {
			// TODO: simplify this
			// Is this even needed? What about default copy constructor?
			p=rhs.p;
			p_center=rhs.p_center;
			T=rhs.T;
			T_center=rhs.T_center;
			rho=rhs.rho;
			a=rhs.a;
			H=rhs.H;
			V=rhs.V;
			V_center=rhs.V_center;
			Vn=rhs.Vn;
			volume=rhs.volume;
			copy(rhs.update.begin(),rhs.update.end(),update.begin());
			return *this;
		}
};

class NS_Face_State {
	public:
		int index;
		double p,T,mu,lambda;
		Vec3D V;
		Vec3D gradu,gradv,gradw,gradT;
		Vec3D normal,tangent1,tangent2,left2right;
		double area;
		int bc;
};

class NS_Fluxes {
	public:
		vector<double> convective,diffusive;
		NS_Fluxes &operator= (const NS_Fluxes & rhs) {
			copy(rhs.convective.begin(),rhs.convective.end(),convective.begin());
			copy(rhs.diffusive.begin(),rhs.diffusive.end(),diffusive.begin());
			return *this;
		}
};

#endif
