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
KSP ksp; // linear solver context
PC pc; // preconditioner context
Vec deltaU,rhs; // solution, residual vectors
Mat impOP; // implicit operator matrix
Vec pseudo_right; // Contribution to rhs due to pseudo time stepping
Mat pseudo_time; // Pseudo time stepping terms

void petsc_init(int argc, char *argv[],double rtol,double abstol,int maxits) {

	// Initialize petsc
	PetscInitialize(&argc,&argv,(char *)0,help);
	//Create nonlinear solver context
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	//KSPGetPC(ksp,&pc);
	//PCSetType(pc,PCILU);
	//PCFactorSetMatOrderingType(pc,MATORDERING_RCM);
	//PCFactorSetReuseOrdering(pc,PETSC_TRUE);
	
	VecCreateMPI(PETSC_COMM_WORLD,grid.cellCount*5,grid.globalCellCount*5,&rhs);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs,&deltaU);
	VecDuplicate(rhs,&pseudo_right);
	VecSet(rhs,0.);
	VecSet(deltaU,0.);

	vector<int> diagonal_nonzeros, off_diagonal_nonzeros;
	int nextCellCount;
	
	// Calculate space necessary for matrix memory allocation
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
	
	
	// Calculate space necessary for pseudo time step terms matrix memory allocation
	diagonal_nonzeros.clear(); off_diagonal_nonzeros.clear();
	for (cit=grid.cell.begin();cit!=grid.cell.end();cit++) {
		for (int i=0;i<5;++i) {
			diagonal_nonzeros.push_back(5);
			off_diagonal_nonzeros.push_back(0);
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
      			&pseudo_time);
	
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetTolerances(ksp,rtol,abstol,1.e15,maxits);
	//KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
	KSPSetType(ksp,KSPFGMRES);
	KSPGMRESSetRestart(ksp,100);
	KSPSetFromOptions(ksp);
	
	return;
} // end petsc_init

void petsc_solve(int &nIter,double &rNorm) {


	
	MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);
	
	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);
	
	if (ps_timeStepMax>1) {
		MatAssemblyBegin(pseudo_time,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(pseudo_time,MAT_FINAL_ASSEMBLY);
		
		VecAssemblyBegin(pseudo_right);
		VecAssemblyEnd(pseudo_right);
		
		VecAXPY(rhs,-1.,pseudo_right); // rhs-=pseudo_right
		MatAXPY(impOP,1.,pseudo_time,DIFFERENT_NONZERO_PATTERN); // impOP+=pseudo_time
	}
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSolve(ksp,rhs,deltaU);
	
	KSPGetIterationNumber(ksp,&nIter);
	KSPGetResidualNorm(ksp,&rNorm); 
	
	if (ps_timeStepMax>1) MatAXPY(impOP,-1.,pseudo_time,DIFFERENT_NONZERO_PATTERN); // impOP-=pseudo_time
	
	int index;
	for (int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<5;++i) {
			index=(grid.myOffset+c)*5+i;
			VecGetValues(deltaU,1,&index,&grid.cell[c].update[i]);
		}
	}

	VecSet(rhs,0.);

	return;
} // end petsc_solve

void petsc_finalize(void) {
	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	PetscFinalize();
	return;
} // end petsc finalize
