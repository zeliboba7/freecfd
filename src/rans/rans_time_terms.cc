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
#include "rans.h"

void RANS::time_terms() {

	PetscInt row,col;
	PetscScalar value;
	
	if (ps_step_max>1) {
		if (ps_step==1) {
			for (int c=0;c<grid[gid].cellCount;++c) {
				row=(grid[gid].myOffset+c)*2;
				value=k.cell(c);
				VecSetValues(soln_n,1,&row,&value,INSERT_VALUES);
				row++; value=omega.cell(c); VecSetValues(soln_n,1,&row,&value,INSERT_VALUES);
			}
			VecAssemblyBegin(soln_n); VecAssemblyEnd(soln_n);
		} else if (ps_step>1) {
			for (int c=0;c<grid[gid].cellCount;++c) {
				row=(grid[gid].myOffset+c)*2;
				value=k.cell(c);
				VecSetValues(pseudo_delta,1,&row,&value,INSERT_VALUES);
				row++; value=omega.cell(c); VecSetValues(pseudo_delta,1,&row,&value,INSERT_VALUES);
			} // pseudo_delta=soln_k
			VecAXPY(pseudo_delta,-1.,soln_n); // pseudo_delta-=soln_n
			VecAssemblyBegin(pseudo_delta); VecAssemblyEnd(pseudo_delta);
			VecSet(pseudo_right,0.);
		}
	}
	
	for (int c=0;c<grid[gid].cellCount;++c) {

		double ps_delta;
		
		// Insert unsteady term
		
		for (int i=0;i<2;++i) {
			row=(grid[gid].myOffset+c)*2+i;
			value=ns[gid].rho.cell(c)*grid[gid].cell[c].volume/dt[gid].cell(c);
			MatSetValues(impOP,1,&row,1,&row,&value,ADD_VALUES);
			if (ps_step>1) {
				VecGetValues(pseudo_delta,1,&row,&ps_delta);
				value*=ps_delta;
				VecSetValues(pseudo_right,1,&row,&value,ADD_VALUES);
			}
		}


		if (ps_step_max>1) {
			for (int i=0;i<2;++i) {
				row=(grid[gid].myOffset+c)*2+i;
				value=ns[gid].rho.cell(c)*grid[gid].cell[c].volume/dtau[gid].cell(c);
				MatSetValues(impOP,1,&row,1,&row,&value,ADD_VALUES);
			}
		}
		
	}
	
	if (ps_step>1) {
		VecAssemblyBegin(pseudo_right); VecAssemblyEnd(pseudo_right);	
		VecAXPY(rhs,-1.,pseudo_right); // rhs-=pseudo_right
	}
	
	
	return;
}

