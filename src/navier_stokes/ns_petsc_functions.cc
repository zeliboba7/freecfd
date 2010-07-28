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
#include "ns.h"

void NavierStokes::petsc_init(void) {
	
	vector<Cell>::iterator cit;
	vector<int>::iterator it;
	
	//Create nonlinear solver context
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	
	VecCreateMPI(PETSC_COMM_WORLD,grid[gid].cellCount*nVars,grid[gid].globalCellCount*nVars,&rhs);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs,&deltaU);
	VecDuplicate(rhs,&pseudo_right);
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
		for (int i=0;i<nVars;++i) {
			diagonal_nonzeros.push_back( (nextCellCount+1)*nVars);
			off_diagonal_nonzeros.push_back( ((*cit).ghosts.size())*nVars);
		}
	}
	
	MatCreateMPIAIJ(
			PETSC_COMM_WORLD,
   			grid[gid].cellCount*nVars,
 			grid[gid].cellCount*nVars,
   			grid[gid].globalCellCount*nVars,
   			grid[gid].globalCellCount*nVars,
   			0,&diagonal_nonzeros[0],
   			0,&off_diagonal_nonzeros[0],
   			&impOP);
	
	
	// Calculate space necessary for pseudo time step terms matrix memory allocation
	diagonal_nonzeros.clear(); off_diagonal_nonzeros.clear();
	for (cit=grid[gid].cell.begin();cit!=grid[gid].cell.end();cit++) {
		for (int i=0;i<nVars;++i) {
			diagonal_nonzeros.push_back(nVars);
			off_diagonal_nonzeros.push_back(0);
		}
	}
	
	MatCreateMPIAIJ(
			PETSC_COMM_WORLD,
   			grid[gid].cellCount*nVars,
   			grid[gid].cellCount*nVars,
   			grid[gid].globalCellCount*nVars,
   			grid[gid].globalCellCount*nVars,
   			0,&diagonal_nonzeros[0],
   			0,&off_diagonal_nonzeros[0],
      			&pseudo_time);
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetTolerances(ksp,rtol,abstol,1.e15,maxits);
	KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
	//KSPSetType(ksp,KSPGMRES);
	KSPSetType(ksp,KSPFGMRES);
	//KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);
	KSPGMRESSetRestart(ksp,100);
	KSPSetFromOptions(ksp);
	
	return;
} 

void NavierStokes::petsc_solve(int &nIter,double &rNorm) {

	MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);
	
	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);
	
	/*
	if (ps_steps_max>1) {
		MatAssemblyBegin(pseudo_time,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(pseudo_time,MAT_FINAL_ASSEMBLY);
		
		VecAssemblyBegin(pseudo_right);
		VecAssemblyEnd(pseudo_right);
		
		VecAXPY(rhs,-1.,pseudo_right); // rhs-=pseudo_right
		MatAXPY(impOP,1.,pseudo_time,DIFFERENT_NONZERO_PATTERN); // impOP+=pseudo_time
	}
	*/
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSolve(ksp,rhs,deltaU);
	
	KSPGetIterationNumber(ksp,&nIter);
	KSPGetResidualNorm(ksp,&rNorm); 
		
	//if (ps_steps_max>1) MatAXPY(impOP,-1.,pseudo_time,DIFFERENT_NONZERO_PATTERN); // impOP-=pseudo_time
	
	int index;
	for (int c=0;c<grid[gid].cellCount;++c) {
		for (int i=0;i<nVars;++i) {
			index=(grid[gid].myOffset+c)*nVars+i;
			VecGetValues(deltaU,1,&index,&update[i].cell(c));
		}
	}

	VecSet(rhs,0.);

	return;
} 

void NavierStokes::petsc_destroy(void) {
	KSPDestroy(ksp);
	MatDestroy(impOP);
	MatDestroy(pseudo_time);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	VecDestroy(pseudo_right);
	PCDestroy(pc);
	return;
} 
