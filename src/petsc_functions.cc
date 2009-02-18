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
#include "petsc_functions.h"
#include "inputs.h"

extern InputFile input;

static char help[] = "Free CFD\n - A free general purpose computational fluid dynamics code";
KSP ksp,kspTurb; // linear solver context
Vec deltaU,deltaUturb,rhs,rhsTurb; // solution, residual vectors
Mat impOP,impOPturb; // implicit operator matrix

void petsc_init(int argc, char *argv[],double rtol,double abstol,int maxits) {

	// Initialize petsc
	PetscInitialize(&argc,&argv,(char *)0,help);
	PC pc,pcTurb; // preconditioner context
	PetscErrorCode ierr;
	//Create nonlinear solver context
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	VecCreateMPI(PETSC_COMM_WORLD,grid.cellCount*5,grid.globalCellCount*5,&rhs);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs,&deltaU);
	VecSet(rhs,0.);
	VecSet(deltaU,0.);

	vector<int> diagonal_nonzeros, off_diagonal_nonzeros;
	int nextCellCount;
	
	// Calculate space necessary for matrix memory alocation
	for (cit=grid.cell.begin();cit!=grid.cell.end();cit++) {
		nextCellCount=0;
		for (it=(*cit).faces.begin();it!=(*cit).faces.end();it++) {
			if (grid.face[*it].bc==INTERNAL) {
				nextCellCount++;
			}
		}
		for (int i=0;i<5;++i) {
			diagonal_nonzeros.push_back( (nextCellCount+1)*5);
			off_diagonal_nonzeros.push_back( ((*cit).ghosts.size())*5);
		}
	}
	
	MatCreateMPIAIJ(
			PETSC_COMM_WORLD,
   			grid.cellCount*5,
 			grid.cellCount*5,
   			grid.globalCellCount*5,
   			grid.globalCellCount*5,
   			0,&diagonal_nonzeros[0],
   			0,&off_diagonal_nonzeros[0],
   			&impOP);
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetTolerances(ksp,rtol,abstol,1.e10,maxits);
	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
	KSPSetType(ksp,KSPFGMRES);
	KSPGMRESSetRestart(ksp,100);
	KSPSetFromOptions(ksp);
	
	if (TURBULENCE_MODEL!=NONE) {
		KSPCreate(PETSC_COMM_WORLD,&kspTurb);
		VecCreateMPI(PETSC_COMM_WORLD,grid.cellCount*2,grid.globalCellCount*2,&rhsTurb);
		VecSetFromOptions(rhsTurb);
		VecDuplicate(rhsTurb,&deltaUturb);
		VecSet(rhsTurb,0.);
		VecSet(deltaUturb,0.);
		
		diagonal_nonzeros.clear(), off_diagonal_nonzeros.clear();
	
		// Calculate space necessary for matrix memory alocation
		for (cit=grid.cell.begin();cit!=grid.cell.end();cit++) {
			nextCellCount=0;
			for (it=(*cit).faces.begin();it!=(*cit).faces.end();it++) {
				if (grid.face[*it].bc==INTERNAL) {
					nextCellCount++;
				}
			}
			for (int i=0;i<2;++i) {
				diagonal_nonzeros.push_back( (nextCellCount+1)*2);
				off_diagonal_nonzeros.push_back( ((*cit).ghosts.size())*2);
			}
		}
	
		MatCreateMPIAIJ(
				PETSC_COMM_WORLD,
				grid.cellCount*2,
				grid.cellCount*2,
				grid.globalCellCount*2,
				grid.globalCellCount*2,
				0,&diagonal_nonzeros[0],
				0,&off_diagonal_nonzeros[0],
				&impOPturb);
	
		KSPSetOperators(kspTurb,impOPturb,impOPturb,SAME_NONZERO_PATTERN);
		KSPSetTolerances(kspTurb,rtol,abstol,1.e10,maxits);
		KSPSetInitialGuessNonzero(kspTurb,PETSC_TRUE);
		KSPSetInitialGuessKnoll(kspTurb,PETSC_TRUE);
		KSPSetType(kspTurb,KSPFGMRES);
		KSPGMRESSetRestart(kspTurb,100);
		KSPSetFromOptions(kspTurb);
	}
	
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
	
	int index;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<5;++i) {
			index=(grid.myOffset+c)*5+i;
			VecGetValues(deltaU,1,&index,&grid.cell[c].update[i]);
		}
	}

	VecSet(rhs,0.);
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);

	return;
} // end petsc_solve

void petsc_solve_turb(int &nIterTurb, double &rNormTurb) {

	MatAssemblyBegin(impOPturb,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOPturb,MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(rhsTurb);
	VecAssemblyEnd(rhsTurb);

	KSPSolve(kspTurb,rhsTurb,deltaUturb);

	KSPGetIterationNumber(kspTurb,&nIterTurb);
	KSPGetResidualNorm(kspTurb,&rNormTurb); 

	int index;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<2;++i) {
			index=(grid.myOffset+c)*2+i;
			VecGetValues(deltaUturb,1,&index,&grid.cell[c].update_turb[i]);
		}
	}

	VecSet(rhsTurb,0.);
	KSPSetOperators(kspTurb,impOPturb,impOPturb,SAME_NONZERO_PATTERN);	

	return;
	
} // end petsc_solve_turb

void petsc_finalize(void) {
	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	KSPDestroy(kspTurb);
	MatDestroy(impOPturb);
	VecDestroy(rhsTurb);
	VecDestroy(deltaUturb);
	PetscFinalize();
	return;
} // end petsc finalize
