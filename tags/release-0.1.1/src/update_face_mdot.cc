/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

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
#include "commons.h"
#include <cmath>
#include "grid.h"
#include "bc.h"
#include "inputs.h"
#include "state_cache.h"
#include "petsc_functions.h"

extern BC bc;
extern InputFile input;

extern void convective_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]);
extern void left_state_update(Cell_State &left,Face_State &face);
extern void right_state_update(Cell_State &left,Cell_State &right,Face_State &face);
extern void face_geom_update(Face_State &face,unsigned int f);

void update_face_mdot(void) {
	
	Cell_State left,right;
	Face_State face;
	Fluxes flux;
	
	flux.diffusive.resize(5);
	flux.convective.resize(5);
	left.update.resize(5);
	right.update.resize(5);

	flux.convective.resize(5);
	
	// Loop through faces
	for (int f=0;f<grid.faceCount;++f) {

		// Populate the state caches
		face_geom_update(face,f);
		left_state_update(left,face);
		right_state_update(left,right,face);

		// Get unperturbed flux values
		convective_face_flux(left,right,face,&flux.convective[0]);
	} // for faces
	return;
} // end function
