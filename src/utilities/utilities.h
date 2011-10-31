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
#ifndef UTILITIES_H
#define UTILITIES_H

#define EPS 1e-10

#include <string>
#include <sstream>
#include <limits>
#include <fstream>
using namespace std;

#include "vec3d.h"

string int2str(int number);

bool fexists(const char *filename);

bool withinBox(Vec3D point, Vec3D corner_1, Vec3D corner_2);
bool withinCylinder(Vec3D point, Vec3D center, double radius, Vec3D axisDirection, double height);
bool withinSphere(Vec3D point, Vec3D center, double radius);

int gelimd(vector<vector<double> > &a,vector<double> &b,vector<double> &x);

#endif
