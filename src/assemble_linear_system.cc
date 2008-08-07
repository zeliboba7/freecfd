/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#include <cmath>
#include "grid.h"
#include "bc.h"
#include "inputs.h"
#include "petsc_functions.h"

extern Grid grid;
extern BC bc;
extern InputFile input;
extern int rank;
extern double Gamma;
extern double Pref;

void inviscid_face_flux(unsigned int f,string order,string limiter,double flux[], int perturb_state=-1, int perturb_var=0, double perturb_value=0.);
		
void assemble_linear_system(void) {

	double flux[5],fluxPlus[5];
	double epsilon;
	unsigned int parent,neighbor,f;
	int row,col;
	PetscScalar value;

	epsilon=sqrt(std::numeric_limits<double>::epsilon());

	string order=input.section["numericalOptions"].strings["order"];
	string limiter=input.section["numericalOptions"].strings["limiter"];
	
	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

		// Get unperturbed flux values
		inviscid_face_flux(f,order,limiter,flux);

		// Fill in residual (rhs vector)
		for (int i=0;i<5;++i) {
			row=grid.cell[parent].globalId*5+i;
			value=-1.*flux[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==-1) { // TODO what if a ghost face??
				row=grid.cell[neighbor].globalId*5+i;
				value*=-1.;
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}
		
		for (int i=0;i<5;++i) {
			// Perturb each left variable
			inviscid_face_flux(f,order,limiter,fluxPlus,0,i,epsilon);
			// Add change of flux (flux Jacobian) to implicit operator
			for (int j=0;j<5;++j) {
				row=grid.cell[parent].globalId*5+j;
				col=grid.cell[parent].globalId*5+i;
				value=(fluxPlus[j]-flux[j])/epsilon;
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				if (grid.face[f].bc==-1) { // TODO what if a ghost face??
					row=grid.cell[neighbor].globalId*5+j;
					value*=-1.;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				}
			} // for j
			if (grid.face[f].bc==-1) {
				// Perturb each right variable
				inviscid_face_flux(f,order,limiter,fluxPlus,1,i,epsilon);
				// Add change of flux (flux Jacobian) to implicit operator
				for (int j=0;j<5;++j) {
					row=grid.cell[neighbor].globalId*5+j;
					col=grid.cell[neighbor].globalId*5+i;
					value=(flux[j]-fluxPlus[j])/epsilon;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					row=grid.cell[parent].globalId*5+j;
					value*=-1.;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				} // for j
			} // if grid.face[f].bc==-1
		} // for i

	} // for faces

	MatAssemblyBegin(impOP,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FLUSH_ASSEMBLY);

	return;
} // end function


