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
#include "inputs.h"

#define CURVILINEAR 1
#define LSQR 2
#define GREENGAUSS 3



using namespace std;

extern vector<Grid> grid;
extern InputFile input;
extern void lsqr_grad_map(int gid,int c);
extern void curvilinear_grad_map(int gid,int c);

void gradient_maps(int gid) {

	int hex_method,prism_method,other_method;

	if (input.section("grid",0).subsection("gradients").get_string("hexmethod")=="curvilinear") hex_method=CURVILINEAR;
	else if (input.section("grid",0).subsection("gradients").get_string("hexmethod")=="lsqr") hex_method=LSQR;
	else if (input.section("grid",0).subsection("gradients").get_string("hexmethod")=="greengauss") hex_method=GREENGAUSS;

	if (input.section("grid",0).subsection("gradients").get_string("prismmethod")=="curvilinear") prism_method=CURVILINEAR;
	else if (input.section("grid",0).subsection("gradients").get_string("prismmethod")=="lsqr") prism_method=LSQR;
	else if (input.section("grid",0).subsection("gradients").get_string("prismmethod")=="greengauss") prism_method=GREENGAUSS;

	if (input.section("grid",0).subsection("gradients").get_string("othermethod")=="curvilinear") other_method=CURVILINEAR;
	else if (input.section("grid",0).subsection("gradients").get_string("othermethod")=="lsqr") other_method=LSQR;
	else if (input.section("grid",0).subsection("gradients").get_string("othermethod")=="greengauss") other_method=GREENGAUSS;

	for (int c=0;c<grid[gid].cellCount;++c) {
		if (grid[gid].cell[c].nodeCount==8) {
			if (hex_method==CURVILINEAR) curvilinear_grad_map(gid,c);
			else if (hex_method==LSQR) lsqr_grad_map(gid,c);
			else if (hex_method==GREENGAUSS) { 
				// if gradMap is not filled, this will be used automatically, so don't need to do anything here 
			}
		} else if (grid[gid].cell[c].nodeCount==6) {
			if (prism_method==CURVILINEAR) curvilinear_grad_map(gid,c);
			else if (prism_method==LSQR) lsqr_grad_map(gid,c);
			else if (prism_method==GREENGAUSS) { 
				// if gradMap is not filled, this will be used automatically, so don't need to do anything here 
			}
		} else {
			if (other_method==CURVILINEAR) curvilinear_grad_map(gid,c);
			else if (other_method==LSQR) lsqr_grad_map(gid,c);
			else if (other_method==GREENGAUSS) { 
				// if gradMap is not filled, this will be used automatically, so don't need to do anything here 
			}
		}
	}
		
	return;
}
