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
#ifndef INTERPOLATE_H
#define INTERPOLATE_H

// Methods
// Weighted Tri-linear
#define WTLI 0
// Inverse distance weighted
#define IDW 1 

#include <iomanip>
using namespace std;

#include "vec3d.h"
#include "utilities.h"
#include "commons.h"

class Interpolate {
private:
	vector<Vec3D>::iterator vit;
	vector<vector<double> > a3,a4;
	vector<double> b3,b4,weights3,weights4;
	
	Vec3D planeNormal,edge1,edge2,edge3,edge4,p_point,centroid;
	Vec3D basis1,basis2,basis3;
	double area,equilateral_area,ave_edge,skewness;
	double tri_weight,tetra_weight,distance;
	double weightSum,tri_weight_sum,tetra_weight_sum;

	bool interpolate_tetra(void);
	bool interpolate_tri(void);
	bool interpolate_line(void);
	bool interpolate_point(void);
	void sort_stencil(bool is_internal);
	
public:
	// Inputs
	int method;
	int max_stencil_size;
	int dimension;
	vector<Vec3D> stencil;
	vector<int> stencil_indices;
	Vec3D point;
	double skewness_tolerance;
	
	// Output (together with updated stencil)
	vector<double> weights;

	
	void init(void);
	void calculate_weights(bool is_internal);
	void flush(void);
	
};

#endif
