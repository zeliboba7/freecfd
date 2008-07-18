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

static char help[] = "Free CFD\n - A free general purpose computational fluid dynamics code";

KSP ksp; // linear solver context
Vec deltaU,rhs,globalUpdate; // solution, residual vectors
Mat impOP; // implicit operator matrix
VecScatter scatterContext;

void petsc_init(int argc, char *argv[]) {
	// Initialize petsc
	PetscInitialize(&argc,&argv,(char *)0,help);
	PC pc; // preconditioner context
	PetscErrorCode ierr;
	//Create nonlinear solver context
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	VecCreateMPI(PETSC_COMM_WORLD,grid.cellCount*5,grid.globalCellCount*5,&rhs);
	VecCreateSeq(PETSC_COMM_SELF,grid.globalCellCount*5,&globalUpdate);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs,&deltaU);
	VecSet(deltaU,0.);
	VecSet(globalUpdate,0.);

	VecScatterCreateToAll(deltaU,&scatterContext,&globalUpdate);

	//PetscScalar *dU,*ff,value;
	
	MatCreateMPIAIJ(PETSC_COMM_WORLD,grid.cellCount*5,grid.cellCount*5,grid.globalCellCount*5,grid.globalCellCount*5,50,PETSC_NULL,0,PETSC_NULL,&impOP);
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetFromOptions(ksp);

	return;
} // end petsc_init

void petsc_solve(int &nIter,double &rNorm) {
	
	MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);

	KSPSolve(ksp,rhs,deltaU);

	KSPGetIterationNumber(ksp,&nIter);
	KSPGetResidualNorm(ksp,&rNorm);

	VecScatterBegin(scatterContext,deltaU,globalUpdate,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(scatterContext,deltaU,globalUpdate,INSERT_VALUES,SCATTER_FORWARD);

	int index;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<5;++i) {
			index=grid.cell[c].globalId*5+i;
			VecGetValues(globalUpdate,1,&index,&grid.cell[c].flux[i]);
		}
	}
	return;
} // end petsc_solve

void petsc_finalize(void) {
	VecScatterDestroy(scatterContext);
	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	PetscFinalize();
	return;
} // end petsc finalize
