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

#ifndef BC_INTERFACE
#define BC_INTERFACE

#include "grid.h"
#include "inputs.h"
#include "commons.h"

extern vector<Grid> grid;
extern InputFile input;
extern vector<int> equations;

class BC_Interface {
public:
	int recv_grid,recv_bc,recv_eqn;
	string recv_var;
	int donor_grid,donor_bc,donor_eqn;
	string donor_var;
	// These two contain essentially a point cloud and the associated data at each point
	vector<double> donor_data;
	vector<Vec3D> donor_point;
	// Current index of the faces in the bc (in order) --> index in donor data
	vector<int> donor_index;
	
	void setup(void);
};

#endif


