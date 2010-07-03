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
#include "utilities.h"

string int2str(int number) {
        std::stringstream ss;
        ss << number;
        return ss.str();
}

bool withinBox(Vec3D point,Vec3D corner_1,Vec3D corner_2) {
	for (int i=0;i<3;++i) {
		if (corner_1[i]>corner_2[i]) {
			if (point[i]>corner_1[i]) return false;
			if (point[i]<corner_2[i]) return false;
		} else {
			if (point[i]<corner_1[i]) return false;
			if (point[i]>corner_2[i]) return false;
		}
	}
	return true;
}

bool withinCylinder(Vec3D point,Vec3D center,double radius,Vec3D axisDirection,double height) {
	
	Vec3D onAxis=(point-center).dot(axisDirection)*axisDirection;
	if (fabs(onAxis)>0.5*height) return false;
	
	Vec3D offAxis=(point-center)-onAxis;
	if (fabs(offAxis)>radius) return false;
	
	return true;
}

bool withinSphere(Vec3D point,Vec3D center,double radius) {
	return (fabs(point-center)<=radius);
}
