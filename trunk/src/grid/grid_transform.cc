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
#include "grid.h"

int Grid::translate(Vec3D begin, Vec3D end) {
	Vec3D diff;
	diff=end-begin;
	for (int n=0;n<globalNodeCount;++n) {
		raw.node[n]+=diff;
	}
	if (Rank==0) cout << "[I grid=" << gid+1 << "] Translated from " << begin << " to " << end << endl;
	return 1;
}

int Grid::scale(Vec3D anchor, double scale) {
	for (int n=0;n<globalNodeCount;++n) {
		raw.node[n]=anchor+scale*(raw.node[n]-anchor);
	}
	if (Rank==0) cout << "[I grid=" << gid+1 << "] scaled by " << scale << " with anchor = " << anchor << endl;
	return 1;
}

int Grid::rotate(Vec3D anchor, Vec3D axis, double angle) {
	// Convert angle to radian
	angle*=4.*atan(1.)/180.;
	// Normalize axis;
	axis=axis.norm();
	Vec3D p;
	for (int n=0;n<globalNodeCount;++n) {
		p=raw.node[n];
		p-=anchor;
		raw.node[n][0]=axis[0]*(axis.dot(p))+(p[0]*(1.-axis[0]*axis[0])-axis[0]*(axis[1]*p[1]+axis[2]*p[2]))*cos(angle)+(-axis[2]*p[1]+axis[1]*p[2])*sin(angle);
		raw.node[n][1]=axis[1]*(axis.dot(p))+(p[1]*(1.-axis[1]*axis[1])-axis[1]*(axis[0]*p[0]+axis[2]*p[2]))*cos(angle)+(axis[2]*p[0]-axis[0]*p[2])*sin(angle);
		raw.node[n][2]=axis[2]*(axis.dot(p))+(p[2]*(1.-axis[2]*axis[2])-axis[2]*(axis[0]*p[0]+axis[1]*p[1]))*cos(angle)+(-axis[1]*p[0]+axis[0]*p[1])*sin(angle);
		raw.node[n]+=anchor;
	}
	angle*=180./(4.*atan(1.));
	if (Rank==0) cout << "[I grid=" << gid+1 << "] rotated by " << angle << " degrees around axis = " << axis << " with anchor = " << anchor << endl;
	return 1;
}
