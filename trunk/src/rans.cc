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
#include "rans.h"

extern double minmod(double a, double b);
extern double doubleMinmod(double a, double b);
extern double harmonic(double a, double b);
extern double superbee(double a, double b);

struct mpiGhost {
	int globalId;
	double vars[3];
};

struct mpiGrad {
	int globalId;
	double grads[6];
};

RANS::RANS(void) {
	kepsilon.sigma_k=1.;
	kepsilon.sigma_omega=0.856;
	kepsilon.beta=0.0828;
	kepsilon.beta_star=0.09;
	kepsilon.kappa=0.41;
	kepsilon.alpha=0.44;
	//kepsilon.alpha=kepsilon.beta/kepsilon.beta_star
	//		-kepsilon.sigma_omega*kepsilon.kappa*kepsilon.kappa/sqrt(kepsilon.beta_star);
	
	komega.sigma_k=0.85;
	komega.sigma_omega=0.5;
	komega.beta=0.075;
	komega.beta_star=0.09;
	komega.kappa=0.41;
	komega.alpha=5./9.;
	//komega.alpha=komega.beta/komega.beta_star
	//		-komega.sigma_omega*komega.kappa*komega.kappa/sqrt(komega.beta_star);
	return; 
} // end RANS::RANS

void RANS::allocate(void) {
	cell.resize(grid.cellCount);
	face.resize(grid.faceCount);
	ghost.resize(grid.ghostCount);
	return; 
} // end RANS::allocate

void RANS::petsc_init(double rtol,double abstol,int maxits) {
	
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
	//KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	KSPSetType(ksp,KSPFGMRES);
	KSPGMRESSetRestart(ksp,100);
	//KSPSetFromOptions(ksp);
	return;
} // end RANS::petsc_init

void RANS::petsc_solve(int &nIter, double &rNorm) {

	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSolve(ksp,rhs,deltaU);

	KSPGetIterationNumber(ksp,&nIter);
	KSPGetResidualNorm(ksp,&rNorm); 

	int index;
	for (int c=0;c<grid.cellCount;++c) {
		for (int i=0;i<2;++i) {
			index=(grid.myOffset+c)*2+i;
			VecGetValues(deltaU,1,&index,&cell[c].update[i]);
		}
	}
	
	return;
	
} // end RANS::petsc_solve

void RANS::petsc_destroy(void) {
	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	return;
} // end RANS::petsc destroy

void RANS::mpi_init(void) {

	// Commit custom communication datatypes
	MPI_Datatype types[2]={MPI_INT,MPI_DOUBLE};
	int block_lengths[2];
	MPI_Aint displacements[2];

	// MPI_GHOST
	mpiGhost dummy1;
	displacements[0]=(long) &dummy1.globalId - (long) &dummy1;
	displacements[1]=(long) &dummy1.vars[0] - (long) &dummy1;
	block_lengths[0]=1;
	block_lengths[1]=3;
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
} // end RANS::mpi_init

void RANS::mpi_update_ghost(void) {
	
	for (int p=0;p<np;++p) {
		if (Rank!=p) {
			mpiGhost sendBuffer[sendCells[p].size()];
			mpiGhost recvBuffer[recvCount[p]];
			int id;
			for (int g=0;g<sendCells[p].size();++g) {
				id=maps.cellGlobal2Local[sendCells[p][g]];
				sendBuffer[g].globalId=grid.cell[id].globalId;
				sendBuffer[g].vars[0]=cell[id].k;
				sendBuffer[g].vars[1]=cell[id].omega;
				sendBuffer[g].vars[2]=cell[id].mu_t;
			}

			MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GHOST,p,0,recvBuffer,recvCount[p],MPI_GHOST,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			for (int g=0;g<recvCount[p];++g) {
				id=maps.ghostGlobal2Local[recvBuffer[g].globalId];
				ghost[id].k=recvBuffer[g].vars[0];
				ghost[id].omega=recvBuffer[g].vars[1];
				ghost[id].mu_t=recvBuffer[g].vars[2];
			}
		}
	}
	
	return;
} // end RANS::mpi_update_ghost

void RANS::mpi_update_ghost_gradients(void) {
	
	// Update ghost gradients of rans variables
	for (int p=0;p<np;++p) {
		mpiGrad sendBuffer[sendCells[p].size()];
		mpiGrad recvBuffer[recvCount[p]];
		int id;
		for (int g=0;g<sendCells[p].size();++g) {
			id=maps.cellGlobal2Local[sendCells[p][g]];
			sendBuffer[g].globalId=grid.cell[id].globalId;
			int count=0;
			for (int var=0;var<2;++var) {
				for (int comp=0;comp<3;++comp) {
					sendBuffer[g].grads[count]=cell[id].grad[var][comp];
					count++;
				}
			}
		}

		MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GRAD,p,0,recvBuffer,recvCount[p],MPI_GRAD,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		for (int g=0;g<recvCount[p];++g) {
			id=maps.ghostGlobal2Local[recvBuffer[g].globalId];
			int count=0;
			for (int var=0;var<2;++var) {
				for (int comp=0;comp<3;++comp) {
					ghost[id].grad[var].comp[comp]=recvBuffer[g].grads[count];
					count++;
				}
			}
		}
	}
	
	return;
} // end RANS::mpi_update_ghost_gradients

void RANS::gradients(void) {

	// Calculate cell gradients
	map<int,Vec3D>::iterator it;
	map<int,double>::iterator fit;
	int f;
	Vec3D areaVec;
	double faceK,faceOmega,faceRho,mu;
	
	for (int c=0;c<grid.cellCount;++c) {
		// Initialize all gradients to zero
		cell[c].grad[0]=0.; cell[c].grad[1]=0.;
		// Add internal and interpartition face contributions
		for (it=grid.cell[c].gradMap.begin();it!=grid.cell[c].gradMap.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				cell[c].grad[0]+=(*it).second*cell[(*it).first].k;
				cell[c].grad[1]+=(*it).second*cell[(*it).first].omega;
			} else { // if contribution is coming from a ghost cell
				cell[c].grad[0]+=(*it).second*ghost[-1*((*it).first+1)].k;
				cell[c].grad[1]+=(*it).second*ghost[-1*((*it).first+1)].omega;
			}
		} // end gradMap loop

 		// Add boundary face contributions
		for (int cf=0;cf<grid.cell[c].faceCount;++cf) {
			f=grid.cell[c].faces[cf];
			if (grid.face[f].bc>=0) { // if a boundary face
				areaVec=grid.face[f].normal*grid.face[f].area/grid.cell[c].volume;
				faceK=0.; faceOmega=0.;
				for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
					if ((*fit).first>=0) { // if contribution is coming from a real cell
						faceK+=(*fit).second*cell[(*fit).first].k;
						faceOmega+=(*fit).second*cell[(*fit).first].omega;
					} else { // if contribution is coming from a ghost cell
						faceK+=(*fit).second*ghost[-1*((*fit).first+1)].k;
						faceOmega+=(*fit).second*ghost[-1*((*fit).first+1)].omega;
					}
				}
				
				faceK=max(faceK,kLowLimit);
				faceK=min(faceK,kHighLimit);
				faceOmega=max(faceOmega,omegaLowLimit);
				
				faceRho=grid.cell[c].rho;
				mu=viscosity;
				
				if (bc.region[grid.face[f].bc].type==INLET) {
					faceK=bc.region[grid.face[f].bc].k;
					faceOmega=bc.region[grid.face[f].bc].omega;
				} else if (bc.region[grid.face[f].bc].type==SYMMETRY) {
					// Symmetry mirrors everything
					faceK=cell[c].k;
					faceOmega=cell[c].omega;
				} else if (bc.region[grid.face[f].bc].type==SLIP) {
					// Slip mirrors everything
					faceK=cell[c].k;
					faceOmega=cell[c].omega;
				} else if (bc.region[grid.face[f].bc].type==NOSLIP) {
					faceK=0.;
					faceOmega=60.*mu/(faceRho*0.075*
							pow(fabs((grid.cell[c].centroid-grid.face[f].centroid).dot(grid.face[f].normal)),2.));
				}
				
				cell[c].grad[0]+=faceK*areaVec;
				cell[c].grad[1]+=faceOmega*areaVec;
				
			} // end if a boundary face
		} // end cell face loop
	} // end cell loop
	
 	
} // end RANS::gradients(void)

void RANS::limit_gradients(void) {
	
	int neighbor,g;
	Vec3D maxGrad[2],minGrad[2];
	if (LIMITER==NONE || order==FIRST) {
		// Don't do anything
	} else {
		for (int c=0;c<grid.cellCount;++c) {
	
			// Initialize min and max to current cells values
			for (int i=0;i<2;++i) maxGrad[i]=minGrad[i]=cell[c].grad[i];
			// Find extremes in neighboring real cells
			for (int cc=0;cc<grid.cell[c].neighborCellCount;++cc) {
				neighbor=grid.cell[c].neighborCells[cc];
				for (int var=0;var<2;++var) {
					for (int comp=0;comp<3;++comp) {
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
			for (int cg=0;cg<grid.cell[c].ghostCount;++cg) {
				g=grid.cell[c].ghosts[cg];
				for (int var=0;var<2;++var) {
					for (int comp=0;comp<3;++comp) {
						maxGrad[var][comp]=max(maxGrad[var][comp]
								,(1.-limiter_sharpening)*ghost[g].grad[var][comp]
								+limiter_sharpening*cell[c].grad[var][comp]);
						minGrad[var][comp]=min(minGrad[var][comp]
								,(1.-limiter_sharpening)*ghost[g].grad[var][comp]
								+limiter_sharpening*cell[c].grad[var][comp]);
					}
				}
			}
			
			if(LIMITER==MINMOD) for (int var=0;var<2;++var) for (int comp=0;comp<3;++comp) cell[c].grad[var][comp]=minmod(maxGrad[var][comp],minGrad[var][comp]);
			if(LIMITER==DOUBLEMINMOD) for (int var=0;var<2;++var) for (int comp=0;comp<3;++comp) cell[c].grad[var][comp]=doubleMinmod(maxGrad[var][comp],minGrad[var][comp]);
			if(LIMITER==HARMONIC) for (int var=0;var<2;++var) for (int comp=0;comp<3;++comp) cell[c].grad[var][comp]=harmonic(maxGrad[var][comp],minGrad[var][comp]);
			if(LIMITER==SUPERBEE) for (int var=0;var<2;++var) for (int comp=0;comp<3;++comp) cell[c].grad[var][comp]=superbee(maxGrad[var][comp],minGrad[var][comp]);
	
		}
	}
	
	return;

} // end RANS::limit_gradients()

void RANS::update(double &resK, double &resOmega) {

	resK=0.; resOmega=0.;
	int counter=0;
	for (int c = 0;c < grid.cellCount;++c) {

		// Limit the update so that k and omega doesn't end up out of limits
		cell[c].update[0]=max(-1.*(cell[c].k-kLowLimit),cell[c].update[0]);
		cell[c].update[0]=min((kHighLimit-cell[c].k),cell[c].update[0]);

		cell[c].update[1]=max(-1.*(cell[c].omega-omegaLowLimit),cell[c].update[1]);
		
		double new_mu_t=grid.cell[c].rho*(cell[c].k+cell[c].update[0])/(cell[c].omega+cell[c].update[1]);
		double mu=viscosity;
		if (new_mu_t/mu>viscosityRatioLimit) {
			counter++; 
			double under_relax;
			double limit_nu=viscosityRatioLimit*mu/grid.cell[c].rho;
			under_relax=(limit_nu*cell[c].omega-cell[c].k)/(cell[c].update[0]-limit_nu*cell[c].update[1]+1.E-8);
			under_relax=0.9*max(1.,under_relax);
			cell[c].update[0]*=under_relax;
			cell[c].update[1]*=under_relax;
		}
		
		cell[c].k += cell[c].update[0];
		cell[c].omega+= cell[c].update[1];
		
		resK+=cell[c].update[0]*cell[c].update[0];
		resOmega+=cell[c].update[1]*cell[c].update[1];

	} // cell loop
	if (counter>0) cout << "[I] Update of k and omega is limited due to viscosityRatioLimit constraint for " << counter << " cells" << endl;
	double residuals[2],totalResiduals[2];
	residuals[0]=resK; residuals[1]=resOmega;
	if (np!=1) {
		MPI_Reduce(&residuals,&totalResiduals,2, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		resK=totalResiduals[0]; resOmega=totalResiduals[1];
	}
	resK=sqrt(resK); resOmega=sqrt(resOmega);
	resK/=resK_norm*double(grid.globalCellCount);
	resOmega/=resOmega_norm*double(grid.globalCellCount);
	
	return;

} // end RANS::update