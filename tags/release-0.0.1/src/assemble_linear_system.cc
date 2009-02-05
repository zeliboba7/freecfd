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
#include "commons.h"
#include <cmath>
#include "grid.h"
#include "bc.h"
#include "inputs.h"
#include "state_cache.h"
#include "petsc_functions.h"

extern BC bc;
extern InputFile input;

bool perturb;

void convective_face_flux(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,double flux[]);
void diffusive_face_flux(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,double flux[]);
void left_state_update(Cell_State &left,Face_State &face,unsigned int f);
void right_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f);
void face_geom_update(Face_State &face,unsigned int f);
void face_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f);
void state_perturb(Cell_State &state,Face_State &face,int var,double epsilon);
void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,int var);

void assemble_linear_system(void) {

	double epsilon;
	unsigned int parent,neighbor,f;
	int row,col;
	PetscScalar value;
	Cell_State left,right;
	Face_State face;

	int nSolVar=5; // Basic equations to solve

	if (TURBULENCE_MODEL!=NONE) nSolVar+=2;

	double flux[nSolVar];
	double conv_flux[nSolVar],conv_fluxPlus[nSolVar];
	double diff_flux[nSolVar],diff_fluxPlus[nSolVar];

	double factor=0.01;

	bool viscousSet=true;
	if (EQUATIONS==EULER) viscousSet=false;
	bool viscous;
	
	bool implicit=true;
	if (TIME_INTEGRATOR==FORWARD_EULER) implicit=false;

	
	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {
		perturb=false;

		viscous=viscousSet;
		
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

		// Populate the state caches
		face_geom_update(face,f);
		left_state_update(left,face,f);
		right_state_update(left,right,face,f);
		if (viscous) face_state_update(left,right,face,f);

		// Flush face fluxes
		for (int m=0;m<7;++m) diff_flux[m]=0.;
		// Get unperturbed flux values
		convective_face_flux(left,right,face,f,conv_flux);
		if (viscous) diffusive_face_flux(left,right,face,f,diff_flux);

		// Integrate boundary fluxes
		if ((timeStep) % integrateBoundaryFreq == 0) {
			if (face.bc>=0) {
				bc.region[face.bc].mass+=conv_flux[0];
				for (int i=0;i<3;++i) bc.region[face.bc].momentum[i]-=(diff_flux[i+1]-conv_flux[i+1]);
				bc.region[face.bc].energy+=conv_flux[4];
			}
		}


		// Fill in rhs vector
		for (int i=0;i<nSolVar;++i) {
			row=(grid.myOffset+parent)*nSolVar+i;
			value=diff_flux[i]-conv_flux[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==INTERNAL) { 
				row=(grid.myOffset+neighbor)*nSolVar+i;
				value*=-1.;
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}

		//  EXPERIMENTAL: blank viscous jacobian entries locally if Cell Reylnolds number is greater than a certain value
// 		if (grid.face[f].bc==-1) {
// 			double cellRe=face.rho*fabs(face.v)*0.5*(grid.cell[parent].lengthScale+grid.cell[neighbor].lengthScale)/mu;
// 			//cout << cellRe << endl;
// 			if (cellRe>10.) {
// 				viscous=false; //cout << "hello" << endl;
// 			}
// 		}

		if (implicit) {
			//if ((timeStep) % jacobianUpdateFreq == 0 || timeStep==restart+1) {
				perturb=true;
		
				// Flux for the diffusive flux jacobian is slightly different
				// due to perturbation, account for that here and
				// recalculate the viscous fluxes
// 				if (viscous) {
// 					for (int m=0;m<7;++m) {
// 						flux[m]=0.;
// 						face_state_adjust(left,right,face,f,m); 
// 					}
// 					diffusive_face_flux(left,right,face,f,diff_flux);
// 				}
		
				for (int i=0;i<nSolVar;++i) {

					if (left.update[i]>0.) {epsilon=max(sqrt_machine_error,factor*left.update[i]);}
					else {epsilon=min(-1.*sqrt_machine_error,factor*left.update[i]); }
					// Perturb left state
					state_perturb(left,face,i,epsilon);
					// If right state is a boundary, correct the condition according to changes in left state
					if (face.bc>=0) right_state_update(left,right,face,f);
					convective_face_flux(left,right,face,f,conv_fluxPlus);
					// Adjust face averages and gradients (crude)
					if (viscous) { // pressure doesn't affect viscous fluxes
						face_state_adjust(left,right,face,f,i);
					 	diffusive_face_flux(left,right,face,f,diff_fluxPlus);
					} else { for (int m=0;m<7;++m) diff_fluxPlus[m]=diff_flux[m]; }
										
					// Restore left state
					state_perturb(left,face,i,-1.*epsilon);
					if (face.bc>=0) right_state_update(left,right,face,f);
					if (viscous) face_state_adjust(left,right,face,f,i);
					// Add change of flux (flux Jacobian) to implicit operator
					for (int j=0;j<nSolVar;++j) {
						row=(grid.myOffset+parent)*nSolVar+j;
						col=(grid.myOffset+parent)*nSolVar+i;
						value=-1.*(diff_fluxPlus[j]-diff_flux[j]-conv_fluxPlus[j]+conv_flux[j])/epsilon;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						if (face.bc==INTERNAL) { 
							row=(grid.myOffset+neighbor)*nSolVar+j;
							value*=-1.;
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						}
					} // for j
		
					if (face.bc==INTERNAL || face.bc==GHOST) { 
						if (right.update[i]>0.) {
							epsilon=max(sqrt_machine_error,factor*right.update[i]); 
						} else {
							epsilon=min(-1.*sqrt_machine_error,factor*right.update[i]); 
						}
						//epsilon=sqrt_machine_error;
						//cout << sqrt_machine_error << "\t" << epsilon << endl;
						
						// Perturb right variable
						state_perturb(right,face,i,epsilon);
						// Adjust face averages and gradients (crude)
						convective_face_flux(left,right,face,f,conv_fluxPlus);
						if (viscous) {
							face_state_adjust(left,right,face,f,i);
							diffusive_face_flux(left,right,face,f,diff_fluxPlus);
						} 
		
						// Restore
						state_perturb(right,face,i,-1.*epsilon);
						if (viscous) face_state_adjust(left,right,face,f,i);
		
						// Add change of flux (flux Jacobian) to implicit operator
						for (int j=0;j<nSolVar;++j) {
							if (face.bc==INTERNAL) { 
								row=(grid.myOffset+neighbor)*nSolVar+j;
								col=(grid.myOffset+neighbor)*nSolVar+i;
								value=(diff_fluxPlus[j]-diff_flux[j]-conv_fluxPlus[j]+conv_flux[j])/epsilon;
								MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
								row=(grid.myOffset+parent)*nSolVar+j;
								value*=-1.;
								MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
							} else  { // Ghost (only add effect on parent cell, effect on itself is taken care of in its own partition
 								row=(grid.myOffset+parent)*nSolVar+j;
 								col=(grid.ghost[-1*neighbor-1].matrix_id)*nSolVar+i;
 								value=(diff_fluxPlus[j]-diff_flux[j]-conv_fluxPlus[j]+conv_flux[j])/epsilon;
 								value*=-1.;
 								MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
							}
						} // for j
					} // if 
				} // for i
			//} // if update Jacobian
		} // if implicit
	} // for faces

// 	if (nSolVar==7) {
// 		for (unsigned int c=0;c<grid.cellCount;++c) {
// 			// Fill in residual (rhs vector)
// 			row=grid.cell[c].globalId*nSolVar+5;
// 			value=-1.*0.09*grid.cell[c].k*grid.cell[c].omega*grid.cell[c].volume;
// 			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
// 			row+=1;
// 			value=-1.*0.3/4.*grid.cell[c].omega*grid.cell[c].omega*grid.cell[c].volume;
// 			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
// 		}
// 	}
	
	return;
} // end function


void left_state_update(Cell_State &left,Face_State &face,unsigned int f) {

	unsigned int parent;
	Vec3D deltaV;
	double delta[7];

	parent=grid.face[f].parent;

	if (order==SECOND) {
		for (unsigned int i=0;i<7;++i) {
			delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].limited_grad[i]);
		}
	} else {
		for (unsigned int i=0;i<7;++i) delta[i]=0.;
	}

	for (unsigned int i=0;i<7;++i) left.update[i]=grid.cell[parent].update[i];
	
	deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];
	
	// Set left primitive variables
	left.p=grid.cell[parent].p+delta[0];
	left.v_center=grid.cell[parent].v;
	left.v=left.v_center+deltaV;
	left.T_center=grid.cell[parent].T;
	left.T=left.T_center+delta[4];
	left.k_center=grid.cell[parent].k;
	left.k=left.k_center+delta[5];
	left.omega_center=grid.cell[parent].omega;
	left.omega=left.omega_center+delta[6];
	left.rho=eos.rho(left.p,left.T);
	left.a=sqrt(Gamma*(left.p+Pref)/left.rho);
	left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
	left.mu=grid.cell[parent].mu;
			
	left.vN[0]=left.v.dot(face.normal);
	left.vN[1]=left.v.dot(face.tangent1);
	left.vN[2]=left.v.dot(face.tangent2);
	
	return;
} // end left_state_update

void right_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f) {

	if (face.bc==INTERNAL) {// internal face

		Vec3D deltaV;
		double delta[7];
		unsigned int neighbor;

  		neighbor=grid.face[f].neighbor;

		if (order==SECOND) {
			for (unsigned int i=0;i<7;++i) {
				delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].limited_grad[i]);
			}
		} else {
			for (unsigned int i=0;i<7;++i) delta[i]=0.;
		}

		for (unsigned int i=0;i<7;++i) right.update[i]=grid.cell[neighbor].update[i];
		
		deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];

		// Set right primitive variables
		right.p=grid.cell[neighbor].p+delta[0];
		right.v_center=grid.cell[neighbor].v;
		right.v=right.v_center+deltaV;
		right.T_center=grid.cell[neighbor].T;
		right.T=right.T_center+delta[4];
		right.k_center=grid.cell[neighbor].k;
		right.k=right.k_center+delta[5];
		right.omega_center=grid.cell[neighbor].omega;
		right.omega=right.omega_center+delta[6];
		right.rho=eos.rho(right.p,right.T);
		right.mu=grid.cell[neighbor].mu;
		
	} else if (face.bc>=0) { // boundary face
		right.mu=left.mu;
		
		for (unsigned int i=0;i<7;++i) right.update[i]=0.;
		
		
		if (bc.region[face.bc].specified==BC_STATE) {
			right.p=bc.region[face.bc].p;
			right.T=bc.region[face.bc].T;
			right.T_center=right.T;
			right.rho=bc.region[face.bc].rho;
		} else if (bc.region[face.bc].specified==BC_P) {
			right.p=0.5*(left.p+bc.region[face.bc].p);
			right.T=left.T; // temperature is extrapolated
			right.T_center=left.T_center+2.*(right.T-left.T_center);
			right.rho=eos.rho(right.p,right.T); 
		} else if (bc.region[face.bc].specified==BC_T) {
			right.T=bc.region[face.bc].T;
			right.T_center=right.T;
			right.p=left.p; // pressure is extrapolated
			right.rho=eos.rho(right.p,right.T); 
		} else if (bc.region[face.bc].specified==BC_RHO) {
			right.rho=0.5*(left.rho+bc.region[face.bc].rho);
			right.p=left.p; // pressure is extrapolated
			right.T=eos.T(right.p,right.rho); 
			right.T_center=left.T_center+2.*(right.T-left.T_center);
		} else {
			// If nothing is specified, everything is extrapolated
			right.p=left.p;
			right.T=left.T; 
			right.T_center=left.T_center+2.*(right.T-left.T_center);
			right.rho=left.rho;
		}
		
		if (bc.region[face.bc].thermalType==ADIABATIC) right.T_center=left.T_center;
		
		if (bc.region[face.bc].type==OUTLET) {
			right.v=left.v;
			right.v_center=left.v_center+2.*(right.v-left.v_center); 
			right.k=left.k;
			right.k_center=left.k_center;
			right.omega=left.omega;
			right.omega_center=left.omega_center;
		} else if (bc.region[face.bc].type==SLIP) {
			right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
			right.v_center=left.v_center-2.*left.v_center.dot(face.normal)*face.normal;
			right.k=left.k;
			right.k_center=left.k_center;
			right.omega=left.omega;
			right.omega_center=left.omega_center;
		} else if (bc.region[face.bc].type==SYMMETRY) {
			right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
			right.v_center=left.v_center-2.*left.v_center.dot(face.normal)*face.normal;
			right.k=left.k;
			right.k_center=left.k_center;
			right.omega=left.omega;
			right.omega_center=left.omega_center;
		} else if (bc.region[face.bc].type==NOSLIP) {
			right.v=-1.*left.v;
			right.v_center=-1.*left.v_center;
			right.k=left.k;
			right.k_center=left.k_center;
			right.omega=left.omega;
			right.omega_center=left.omega_center;
		} else if (bc.region[face.bc].type==INLET) {
			right.v=bc.region[face.bc].v;
			right.v_center=right.v;
			right.k=bc.region[face.bc].k;
			right.omega=bc.region[face.bc].omega;
			right.k_center=right.k;
			right.omega_center=right.omega;
		}

	} else { // partition boundary

		int g=-1*grid.face[f].neighbor-1; // ghost cell index
		Vec3D deltaV;
		double delta[7];
		if (order==SECOND) {
			for (unsigned int i=0;i<7;++i) {
				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].limited_grad[i]);
			}
		}

		for (unsigned int i=0;i<7;++i) right.update[i]=grid.ghost[g].update[i];
		
		deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];
		right.p=grid.ghost[g].p+delta[0];	
		right.v_center=grid.ghost[g].v;
		right.v=right.v_center+deltaV;
		right.T_center=grid.ghost[g].T;
		right.T=right.T_center+delta[4];
		right.k=grid.ghost[g].k+delta[5];
		right.omega=grid.ghost[g].omega+delta[6];
		right.rho=eos.rho(right.p,right.T);
		right.mu=grid.ghost[g].mu;
	}

	right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
	right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
	right.vN[0]=right.v.dot(face.normal);
	right.vN[1]=right.v.dot(face.tangent1);
	right.vN[2]=right.v.dot(face.tangent2);
	
	return;
} // end right_state_update

void face_geom_update(Face_State &face,unsigned int f) {
	face.normal=grid.face[f].normal;
	face.tangent1=(grid.face[f].centroid-grid.face[f].node(0));
	face.tangent1/=fabs(face.tangent1);
	// Cross the tangent vector with the normal vector to get the second tangent
	face.tangent2=(face.normal).cross(face.tangent1);
	face.area=grid.face[f].area;
	face.bc=grid.face[f].bc;
	Vec3D leftCentroid,rightCentroid;
	leftCentroid=grid.cell[grid.face[f].parent].centroid;
	if (face.bc==INTERNAL) { rightCentroid=grid.cell[grid.face[f].neighbor].centroid;}
	else if (face.bc>=0) {
		rightCentroid=leftCentroid+2.*(grid.face[f].centroid-leftCentroid).dot(face.normal)*face.normal;
	} else { rightCentroid=grid.ghost[-1*grid.face[f].neighbor-1].centroid;}
	face.left2right=rightCentroid-leftCentroid;
	return;
} // end face_geom_update

void face_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f) {

	map<int,double>::iterator fit;
	// Find face averaged variables
	face.gradU=0.; face.gradV=0.; face.gradW=0.; face.gradT=0.; face.gradK=0.; face.gradOmega=0.;
	for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			face.gradU+=(*fit).second*grid.cell[(*fit).first].grad[1];
			face.gradV+=(*fit).second*grid.cell[(*fit).first].grad[2];
			face.gradW+=(*fit).second*grid.cell[(*fit).first].grad[3];
			face.gradT+=(*fit).second*grid.cell[(*fit).first].grad[4];
			face.gradK+=(*fit).second*grid.cell[(*fit).first].grad[5];
			face.gradOmega+=(*fit).second*grid.cell[(*fit).first].grad[6];
		} else { // if contribution is coming from a ghost cell
			face.gradU+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[1];
			face.gradV+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[2];
			face.gradW+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[3];
			face.gradT+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[4];
			face.gradK+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[5];
			face.gradOmega+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[6];
		}
	}
	
		
	face.gradU-=face.gradU.dot(face.normal)*face.normal;
	face.gradU+=((right.v_center[0]-left.v_center[0])/(face.left2right.dot(face.normal)))*face.normal;

	face.gradV-=face.gradV.dot(face.normal)*face.normal;
	face.gradV+=((right.v_center[1]-left.v_center[1])/(face.left2right.dot(face.normal)))*face.normal;

	face.gradW-=face.gradW.dot(face.normal)*face.normal;
	face.gradW+=((right.v_center[2]-left.v_center[2])/(face.left2right.dot(face.normal)))*face.normal;
	//cout << face.left2right.norm() << endl;
	// Boundary conditions are already taken care of in right state update
	face.p=0.5*(left.p+right.p);
	face.v=0.5*(left.v+right.v);
	face.T=0.5*(left.T+right.T);
	face.k=0.5*(left.k+right.k);
	face.omega=0.5*(left.omega+right.omega);
	face.rho=0.5*(left.rho+right.rho);
	face.mu=0.5*(left.mu+right.mu);
	
	return;
} // end face_state_update

void state_perturb(Cell_State &state,Face_State &face,int var,double epsilon) {

	switch (var)
	{
		case 0 : // p
			state.p+=epsilon;
			state.rho=eos.rho(state.p,state.T);
			state.a=sqrt(Gamma*(state.p+Pref)/state.rho);
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			break;
		case 1 : // u
			state.v[0]+=epsilon;
			state.v_center[0]+=epsilon;
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 2 : // v
			state.v[1]+=epsilon;
			state.v_center[1]+=epsilon;
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 3 : // w
			state.v[2]+=epsilon;
			state.v_center[2]+=epsilon;
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 4 : // T
			state.T+=epsilon;
			state.T_center+=epsilon;
			state.rho=eos.rho(state.p,state.T);
			state.a=sqrt(Gamma*(state.p+Pref)/state.rho);
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			break;
		case 5 : // k
			state.k+=epsilon;
			state.k_center+=epsilon;
			break;
		case 6 : // omega
			state.omega+=epsilon;
			state.omega_center+=epsilon;
			break;
	}
	
	return;
} // end state_perturb

void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,int var) {

	double distance;
	Vec3D direction;
	
	if (var!=0) {
		distance=fabs(face.left2right);
		direction=face.left2right/distance;
	}

	switch (var)
	{
		case 0 : // p
			face.p=0.5*(left.p+right.p);
			face.rho=eos.rho(face.p,face.T);
			break;
		case 1 : // u
			face.v[0]=0.5*(left.v[0]+right.v[0]);
// 			face.gradU-=face.gradU.dot(direction)*direction;
// 			face.gradU+=((right.v_center[0]-left.v_center[0])/distance)*direction;
			face.gradU-=face.gradU.dot(face.normal)*face.normal;
			face.gradU+=((right.v_center[0]-left.v_center[0])/(face.left2right.dot(face.normal)))*face.normal;

			break;
		case 2 : // v
			face.v[1]=0.5*(left.v[1]+right.v[1]);
			face.gradV-=face.gradV.dot(face.normal)*face.normal;
			face.gradV+=((right.v_center[1]-left.v_center[1])/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 3 : // w
			face.v[2]=0.5*(left.v[2]+right.v[2]);
			face.gradW-=face.gradW.dot(face.normal)*face.normal;
			face.gradW+=((right.v_center[2]-left.v_center[2])/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 4 : // T
			face.T=0.5*(left.T+right.T);
			face.rho=eos.rho(face.p,face.T);
			face.gradT-=face.gradT.dot(direction)*direction;
			face.gradT+=((right.T_center-left.T_center)/distance)*direction;
			break;
		case 5 : // k
			face.k=0.5*(left.k+right.k);
			face.gradK-=face.gradK.dot(direction)*direction;
			face.gradK+=((right.k_center-left.k_center)/distance)*direction;
			break;
		case 6 : // omega
			face.omega=0.5*(left.omega+right.omega);
			face.gradOmega-=face.gradOmega.dot(direction)*direction;
			face.gradOmega+=((right.omega_center-left.omega_center)/distance)*direction;
			break;
	}

	return;
} // end face_state_adjust
