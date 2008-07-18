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
#include "petsc_functions.h"

extern double Gamma,dt;

void linear_system_initialize() {

	MatZeroEntries(impOP);

	PetscInt counter=0;

	// Variables for preconditioner
	double Beta,Mach,Mach2,soundSpeed2;
	double P51,P52,P53,P54,P55,d;
	PetscInt row,col,index;
	PetscScalar value;
			
	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<5;++i) {
			index=grid.cell[c].globalId*5+i;

			// Fill in right hand side vector
			value=-1.*grid.cell[c].flux[i];
			VecSetValues(rhs,1,&index,&value,INSERT_VALUES);

			// Preconditioner
			soundSpeed2=Gamma*grid.cell[c].p/grid.cell[c].rho;
			Mach2=grid.cell[c].v.dot(grid.cell[c].v)/soundSpeed2;
			Mach=sqrt(Mach2);
			if (Mach<=1.e-5) {
				Mach=1.e-5;
			} else if (Mach<1.) {
			// use local
			} else {
				Mach=1.;
			}

			Mach2=Mach*Mach;
			P51=grid.cell[c].v.dot(grid.cell[c].v)/3.*(Mach2-1.);
			P52=grid.cell[c].v.comp[0]*(1.-1./Mach2);
			P53=grid.cell[c].v.comp[1]*(1.-1./Mach2);
			P54=grid.cell[c].v.comp[2]*(1.-1./Mach2);
			P55=1./Mach2;

			d=grid.cell[c].volume/dt;
			
			if (i==5) {
				value=P51*d;
				col=index-4;
				MatSetValues(impOP,1,&index,1,&col,&value,INSERT_VALUES);
				value=P52*d;
				col=index-3;
				MatSetValues(impOP,1,&index,1,&col,&value,INSERT_VALUES);
				value=P53*d;
				col=index-2;
				MatSetValues(impOP,1,&index,1,&col,&value,INSERT_VALUES);
				value=P54*d;
				col=index-1;
				MatSetValues(impOP,1,&index,1,&col,&value,INSERT_VALUES);
				value=P55*d;
				col=index;
				MatSetValues(impOP,1,&index,1,&col,&value,INSERT_VALUES);
			} else {
				value=grid.cell[c].volume/dt;
				MatSetValues(impOP,1,&index,1,&index,&value,INSERT_VALUES);
			}
		}
	}

	MatAssemblyBegin(impOP,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FLUSH_ASSEMBLY);

	
	return;
}
