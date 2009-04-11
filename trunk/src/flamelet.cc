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
#include "petscksp.h"
#include "flamelet.h"
#include "rans.h"

extern double minmod(double a, double b);
extern double doubleMinmod(double a, double b);
extern double harmonic(double a, double b);
extern double superbee(double a, double b);
extern RANS rans;

struct mpiGhost {
	unsigned int globalId;
	double vars[7]; // Z, Zvar,Chi,mu,diffusivity,R and gamma
};

struct mpiGrad {
	unsigned int globalId;
	double grads[6];
};

Flamelet::Flamelet(void) {
	constants.sigma_t=0.7;
	constants.Cg=2.86;
	constants.Cd=2.0;
	return; 
} // end Flamelet::Flamelet

void Flamelet::allocate(void) {
	cell.resize(grid.cellCount);
	face.resize(grid.faceCount);
	ghost.resize(grid.ghostCount);
	return; 
} // end Flamelet::allocate

void Flamelet::petsc_init(double rtol,double abstol,int maxits) {
	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	VecCreateMPI(PETSC_COMM_WORLD,grid.cellCount*2,grid.globalCellCount*2,&rhs);
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
    		&impOP);
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetTolerances(ksp,rtol,abstol,1.e10,maxits);
	KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
	KSPSetType(ksp,KSPFGMRES);
	KSPGMRESSetRestart(ksp,100);
	KSPSetFromOptions(ksp);
	return;
} // end Flamelet::petsc_init

void Flamelet::petsc_solve(int &nIter, double &rNorm) {

	MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);

	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);	
	KSPSolve(ksp,rhs,deltaU);

	KSPGetIterationNumber(ksp,&nIter);
	KSPGetResidualNorm(ksp,&rNorm); 

	int index;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<2;++i) {
			index=(grid.myOffset+c)*2+i;
			VecGetValues(deltaU,1,&index,&cell[c].update[i]);
		}
		//cell[c].update[0]=0.;
	}

	return;
	
} // end Flamelet::petsc_solve

void Flamelet::petsc_destroy(void) {
	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	return;
} // end Flamelet::petsc destroy

void Flamelet::mpi_init(void) {

	// Commit custom communication datatypes
	MPI_Datatype types[2]={MPI_UNSIGNED,MPI_DOUBLE};
	int block_lengths[2];
	MPI_Aint displacements[2];

	// MPI_GHOST
	mpiGhost dummy1;
	displacements[0]=(long) &dummy1.globalId - (long) &dummy1;
	displacements[1]=(long) &dummy1.vars[0] - (long) &dummy1;
	block_lengths[0]=1;
	block_lengths[1]=7;
	MPI_Type_create_struct(2,block_lengths,displacements,types,&MPI_GHOST);
	MPI_Type_commit(&MPI_GHOST);

	// MPI_GRAD
	mpiGrad dummy2;
	displacements[0]=(long) &dummy2.globalId - (long) &dummy2;
	displacements[1]=(long) &dummy2.grads[0] - (long) &dummy2;
	block_lengths[0]=1;
	block_lengths[1]=6;
	MPI_Type_create_struct(2,block_lengths,displacements,types,&MPI_GRAD);
	MPI_Type_commit(&MPI_GRAD);

	return;
} // end Flamelet::mpi_init

void Flamelet::mpi_update_ghost(void) {
	
	for (unsigned int p=0;p<np;++p) {
		if (Rank!=p) {
			mpiGhost sendBuffer[sendCells[p].size()];
			mpiGhost recvBuffer[recvCount[p]];
			int id;
			for (unsigned int g=0;g<sendCells[p].size();++g) {
				id=maps.cellGlobal2Local[sendCells[p][g]];
				sendBuffer[g].globalId=grid.cell[id].globalId;
				sendBuffer[g].vars[0]=cell[id].Z;
				sendBuffer[g].vars[1]=cell[id].Zvar;
				sendBuffer[g].vars[2]=cell[id].Chi;
				sendBuffer[g].vars[3]=cell[id].mu;
				sendBuffer[g].vars[4]=cell[id].diffusivity;
				sendBuffer[g].vars[5]=cell[id].R;
				sendBuffer[g].vars[6]=cell[id].gamma;
			}

			MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GHOST,p,0,recvBuffer,recvCount[p],MPI_GHOST,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			for (unsigned int g=0;g<recvCount[p];++g) {
				id=maps.ghostGlobal2Local[recvBuffer[g].globalId];
				ghost[id].Z=recvBuffer[g].vars[0];
				ghost[id].Zvar=recvBuffer[g].vars[1];
				ghost[id].Chi=recvBuffer[g].vars[2];
				ghost[id].mu=recvBuffer[g].vars[3];
				ghost[id].diffusivity=recvBuffer[g].vars[4];
				ghost[id].R=recvBuffer[g].vars[5];
				ghost[id].gamma=recvBuffer[g].vars[6];
			}
		}
	}
	
	return;
} // end Flamelet::mpi_update_ghost

void Flamelet::mpi_update_ghost_gradients(void) {
	
	for (unsigned int p=0;p<np;++p) {
		mpiGrad sendBuffer[sendCells[p].size()];
		mpiGrad recvBuffer[recvCount[p]];
		int id;
		for (unsigned int g=0;g<sendCells[p].size();++g) {
			id=maps.cellGlobal2Local[sendCells[p][g]];
			sendBuffer[g].globalId=grid.cell[id].globalId;
			int count=0;
			for (unsigned int var=0;var<2;++var) {
				for (unsigned int comp=0;comp<3;++comp) {
					sendBuffer[g].grads[count]=cell[id].grad[var][comp];
					count++;
				}
			}
		}

		MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GRAD,p,0,recvBuffer,recvCount[p],MPI_GRAD,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		for (unsigned int g=0;g<recvCount[p];++g) {
			id=maps.ghostGlobal2Local[recvBuffer[g].globalId];
			int count=0;
			for (unsigned int var=0;var<2;++var) {
				for (unsigned int comp=0;comp<3;++comp) {
					ghost[id].grad[var][comp]=recvBuffer[g].grads[count];
					count++;
				}
			}
		}
	}
	
	return;
} // end Flamelet::mpi_update_ghost_gradients

void Flamelet::gradients(void) {

	// Calculate cell gradients
	map<int,Vec3D>::iterator it;
	map<int,double>::iterator fit;
	unsigned int f;
	Vec3D areaVec;
	double faceZ,faceZvar;
	
	for (unsigned int c=0;c<grid.cellCount;++c) {
		// Initialize all gradients to zero
		cell[c].grad[0]=0.; cell[c].grad[1]=0.;
		// Add internal and interpartition face contributions
		for (it=grid.cell[c].gradMap.begin();it!=grid.cell[c].gradMap.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				cell[c].grad[0]+=(*it).second*cell[(*it).first].Z;
				cell[c].grad[1]+=(*it).second*cell[(*it).first].Zvar;
			} else { // if contribution is coming from a ghost cell
				cell[c].grad[0]+=(*it).second*ghost[-1*((*it).first+1)].Z;
				cell[c].grad[1]+=(*it).second*ghost[-1*((*it).first+1)].Zvar;
			}
		} // end gradMap loop

 		// Add boundary face contributions
		for (unsigned int cf=0;cf<grid.cell[c].faceCount;++cf) {
			f=grid.cell[c].faces[cf];
			if (grid.face[f].bc>=0) { // if a boundary face
				areaVec=grid.face[f].normal*grid.face[f].area/grid.cell[c].volume;
				
				faceZ=cell[c].Z;
				faceZvar=cell[c].Zvar;
				if (bc.region[grid.face[f].bc].type==INLET) {
					faceZ=bc.region[grid.face[f].bc].Z;
					faceZvar=bc.region[grid.face[f].bc].Zvar;
				}
				
				cell[c].grad[0]+=faceZ*areaVec;
				cell[c].grad[1]+=faceZvar*areaVec;
				
			} // end if a boundary face
		} // end cell face loop
	} // end cell loop
	
 	
} // end Flamelet::gradients(void)

void Flamelet::limit_gradients(void) {
	
	unsigned int neighbor,g;
	Vec3D maxGrad[2],minGrad[2];
	if (LIMITER==NONE || order==FIRST) {
		// Don't do anything
	} else {
		for (unsigned int c=0;c<grid.cellCount;++c) {
	
			// Initialize min and max to current cells values
			for (unsigned int i=0;i<2;++i) maxGrad[i]=minGrad[i]=cell[c].grad[i];
			// Find extremes in neighboring real cells
			for (unsigned int cc=0;cc<grid.cell[c].neighborCellCount;++cc) {
				neighbor=grid.cell[c].neighborCells[cc];
				for (unsigned int var=0;var<2;++var) {
					for (unsigned int comp=0;comp<3;++comp) {
						maxGrad[var][comp]=max(maxGrad[var][comp]
								,(1.-limiter_sharpening)*cell[neighbor].grad[var][comp]
								+limiter_sharpening*cell[c].grad[var][comp]);
						
						minGrad[var][comp]=min(minGrad[var][comp]
								,(1.-limiter_sharpening)*cell[neighbor].grad[var][comp]
								+limiter_sharpening*cell[c].grad[var][comp]);
					}
				}
			}
			// Find extremes in neighboring ghost cells
			for (unsigned int cg=0;cg<grid.cell[c].ghostCount;++cg) {
				g=grid.cell[c].ghosts[cg];
				for (unsigned int var=0;var<2;++var) {
					for (unsigned int comp=0;comp<3;++comp) {
						maxGrad[var][comp]=max(maxGrad[var][comp]
								,(1.-limiter_sharpening)*ghost[g].grad[var][comp]
								+limiter_sharpening*cell[c].grad[var][comp]);
						minGrad[var][comp]=min(minGrad[var][comp]
								,(1.-limiter_sharpening)*ghost[g].grad[var][comp]
								+limiter_sharpening*cell[c].grad[var][comp]);
					}
				}
			}
			
			if(LIMITER==MINMOD) for (unsigned int var=0;var<2;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].grad[var][comp]=minmod(maxGrad[var][comp],minGrad[var][comp]);
			if(LIMITER==DOUBLEMINMOD) for (unsigned int var=0;var<2;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].grad[var][comp]=doubleMinmod(maxGrad[var][comp],minGrad[var][comp]);
			if(LIMITER==HARMONIC) for (unsigned int var=0;var<2;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].grad[var][comp]=harmonic(maxGrad[var][comp],minGrad[var][comp]);
			if(LIMITER==SUPERBEE) for (unsigned int var=0;var<2;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].grad[var][comp]=superbee(maxGrad[var][comp],minGrad[var][comp]);
	
		}
	}
	
	return;

} // end Flamelet::limit_gradients()

void Flamelet::update(double &resZ, double &resZvar) {
	

	resZ=0.; resZvar=0.;
	
	for (unsigned int c = 0;c < grid.cellCount;++c) {
		
		// Limit the update so Z doesn't end up negative
		cell[c].update[0]=max(-1.*cell[c].Z,cell[c].update[0]);
		// Limit the update so Zvar doesn't end up negative
		cell[c].update[1]=max(-1.*cell[c].Zvar,cell[c].update[1]);
		// Limit the update so Z doesn't end up larger than 1
		cell[c].update[0]=min(1.-cell[c].Z,cell[c].update[0]);
		// Limit the update so Z doesn't end up larger than 0.25
		cell[c].update[1]=min(0.25-cell[c].Zvar,cell[c].update[1]);
		
		//cell[c].Z+=cell[c].update[0];
		cell[c].Zvar+=cell[c].update[1];
		
		//resZ+=cell[c].update[0]*cell[c].update[0];
		resZvar+=cell[c].update[1]*cell[c].update[1];
		
	} // cell loop

	double residuals[2],totalResiduals[2];
	residuals[0]=resZ; residuals[1]=resZvar;
	if (np!=1) {
		MPI_Reduce(&residuals,&totalResiduals,2, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		//resZ=totalResiduals[0]; 
		resZvar=totalResiduals[1];
	}
	//resZ=sqrt(resZ); 
	resZvar=sqrt(resZvar);
	
	return;

} // end Flamelet::update

void Flamelet::lookup(void) {
	
	double new_T;
	
	for (unsigned int c=0;c<grid.cellCount;++c) {
		
		cell[c].Chi=log10(2.0*rans.kepsilon.beta_star*rans.cell[c].omega*cell[c].Zvar);	
		table.get_rho_T_comp(grid.cell[c].p,cell[c].Z,cell[c].Zvar,cell[c].Chi,grid.cell[c].rho,new_T);
		grid.cell[c].update[4]=new_T-grid.cell[c].T;
		grid.cell[c].T=new_T;
		
		cell[c].R=UNIV_GAS_CONST/table.get_Mw(cell[c].Z,cell[c].Zvar,cell[c].Chi);
		cell[c].gamma=table.get_gamma(cell[c].Z,cell[c].Zvar,cell[c].Chi,false);
		cell[c].mu=table.get_mu(cell[c].Z,cell[c].Zvar,cell[c].Chi,false);
		cell[c].diffusivity=table.get_diffusivity(cell[c].Z,cell[c].Zvar,cell[c].Chi,false);
		
	} // cell loop
	
	return;

} // end Flamelet::lookup

void Flamelet::update_face_properties(void) {

	for (unsigned int f=0;f<grid.faceCount;++f) {

		unsigned int parent,neighbor;
		double leftZ,leftZvar,rightZ,rightZvar;
		Vec3D faceGradZ, faceGradZvar, left2right;
		bool extrapolated;

		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		
		// Get left, right and face values of Z and Zvar as well as the face normal gradients
		get_Z_Zvar(parent,neighbor,f,
			   leftZ,leftZvar,rightZ,rightZvar,
      			   faceGradZ,faceGradZvar,left2right,extrapolated);
		
		double Chi=cell[parent].Chi;
		if (grid.face[f].bc==INTERNAL) {
			Chi=0.5*(Chi+cell[neighbor].Chi);
		} else if (grid.face[f].bc==GHOST) {
			Chi=0.5*(Chi+ghost[-1*neighbor-1].Chi);
		}
		
		face[f].RL=UNIV_GAS_CONST/table.get_Mw(leftZ,leftZvar,Chi);
		face[f].gammaL=table.get_gamma(leftZ,leftZvar,Chi,false);
		face[f].mu=table.get_mu(leftZ,leftZvar,Chi,false);
		face[f].RR=UNIV_GAS_CONST/table.get_Mw(rightZ,rightZvar,Chi);
		face[f].gammaR=table.get_gamma(rightZ,rightZvar,Chi,false);
		face[f].mu+=table.get_mu(rightZ,rightZvar,Chi,false);
		face[f].mu*=0.5;
	}
	
	return;
}

