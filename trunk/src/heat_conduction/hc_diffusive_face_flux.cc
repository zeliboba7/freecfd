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
#include "hc.h"

void HeatConduction::diffusive_face_flux(HC_Face_State &face,double &flux) {

	Vec3D areaVec;

	areaVec=face.normal*face.area;
	
	if (face.bc>=0) {
		if (bc[gid][face.bc].thermalType==FIXED_Q) {
			flux=qdot.bc(face.bc,face.index)*face.area;
		} else if (bc[gid][face.bc].thermalType==FIXED_T) {
			flux=face.lambda/(material.density*material.Cp(face.T))*face.gradT.dot(areaVec);
		} else {
			// Adiabatic
			flux=0.;
		}
	} else {
		flux=face.lambda/(material.density*material.Cp(face.T))*face.gradT.dot(areaVec);
	}
	
	return;
}

