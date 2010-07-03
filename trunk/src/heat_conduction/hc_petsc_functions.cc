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
#include "hc.h"

void HeatConduction::petsc_init(void) {
	
	vector<Cell>::iterator cit;
	vector<int>::iterator it;
	
	//Create nonlinear solver context
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	
	VecCreateMPI(PETSC_COMM_WORLD,grid[gid].cellCount,grid[gid].globalCellCount,&rhs);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs,&deltaU);
	VecSet(rhs,0.);
	VecSet(deltaU,0.);

	vector<int> diagonal_nonzeros, off_diagonal_nonzeros;
	int nextCellCount;
	
	// Calculate space necessary for matrix memory allocation
	for (cit=grid[gid].cell.begin();cit!=grid[gid].cell.end();cit++) {
		nextCellCount=0;
		for (it=(*cit).faces.begin();it!=(*cit).faces.end();it++) {
			if (grid[gid].face[*it].bc==INTERNAL_FACE) {
				nextCellCount++;
			}
		}
		diagonal_nonzeros.push_back((nextCellCount+1));
		off_diagonal_nonzeros.push_back(((*cit).ghosts.size()));
	}
	
	MatCreateMPIAIJ(
			PETSC_COMM_WORLD,
   			grid[gid].cellCount,
 			grid[gid].cellCount,
   			grid[gid].globalCellCount,
   			grid[gid].globalCellCount,
   			0,&diagonal_nonzeros[0],
   			0,&off_diagonal_nonzeros[0],
   			&impOP);

	diagonal_nonzeros.clear(); off_diagonal_nonzeros.clear();
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetTolerances(ksp,rtol,abstol,1.e15,maxits);
	KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
	KSPSetType(ksp,KSPFGMRES);
	KSPGMRESSetRestart(ksp,100);
	KSPSetFromOptions(ksp);
	
	return;
} 

void HeatConduction::petsc_solve(int &nIter,double &rNorm) {

	MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);
	
	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSolve(ksp,rhs,deltaU);
	
	KSPGetIterationNumber(ksp,&nIter);
	KSPGetResidualNorm(ksp,&rNorm); 
		
	int index;
	for (int c=0;c<grid[gid].cellCount;++c) {
		index=grid[gid].myOffset+c;
		VecGetValues(deltaU,1,&index,&update.cell(c));
	}

	VecSet(rhs,0.);

	return;
} 

void HeatConduction::petsc_destroy(void) {
	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	return;
} 
