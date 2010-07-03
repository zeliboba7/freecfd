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
#ifndef HC_STATE_CACHE_H
#define HC_STATE_CACHE_H

class HC_Cell_State {
	public:
		double T,T_center,volume;
		double update;
		HC_Cell_State &operator= (const HC_Cell_State & rhs) {
			// TODO: simplify this
			// Is this even needed? What about default copy constructor?
			T=rhs.T;
			T_center=rhs.T_center;
			volume=rhs.volume;
			update=rhs.update;
			return *this;
		}
};

class HC_Face_State {
	public:
		int index;
		double T,lambda;
		Vec3D gradT;
		Vec3D normal,tangent1,tangent2,left2right;
		double area;
		int bc;
};

#endif
