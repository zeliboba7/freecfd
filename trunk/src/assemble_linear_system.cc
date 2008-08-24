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
#include <cmath>
#include "grid.h"
#include "bc.h"
#include "inputs.h"
#include "state_cache.h"
#include "petsc_functions.h"

extern Grid grid;
extern BC bc;
extern InputFile input;
extern int Rank;
extern double Gamma;
extern double Pref;

string order;
string limiter;

void convective_face_flux(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,double flux[]);
void diffusive_face_flux(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,double flux[]);
void left_state_update(Cell_State &left,Face_State &face,unsigned int f);
void right_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f);
void face_geom_update(Face_State &face,unsigned int f);
void face_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f);
void left_state_perturb(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,int var,double epsilon);
void right_state_perturb(Cell_State &right,Face_State &face,int var,double epsilon);
void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,int var);

void assemble_linear_system(void) {

	double epsilon;
	unsigned int parent,neighbor,f;
	int row,col;
	PetscScalar value;
	Cell_State left,right;
	Face_State face;
	
	int nSolVar=5; // Basic equations to solve

	if (input.section["turbulence"].strings["model"]!="none") nSolVar+=2;

	double flux[nSolVar],fluxPlus[nSolVar];

	epsilon=sqrt(std::numeric_limits<double>::epsilon()); // A small number

	order=input.section["numericalOptions"].strings["order"];
	limiter=input.section["numericalOptions"].strings["limiter"];

	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {

		// TODO Some Ideas:
		// If a slip, symmetry or no-slip face, only look at pressure perturbation in convective fluxes
		// rho and p perturbations doesn't have any effect in diffusive fluxes
		
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

		// Populate the state caches
		face_geom_update(face,f);
		left_state_update(left,face,f);
		right_state_update(left,right,face,f);
		face_state_update(left,right,face,f);

		// Flush face fluxes
		for (int m=0;m<7;++m) flux[m]=0.;
		// Get unperturbed flux values
		convective_face_flux(left,right,face,f,flux);
		diffusive_face_flux(left,right,face,f,flux);

		// Integrate boundary momentum values
		if (face.bc>=0) {
			Vec3D momentum;
			momentum.comp[0]=flux[1];//-(Pref*face.normal.comp[0])*face.area;
			momentum.comp[1]=flux[2];//-(Pref*face.normal.comp[1])*face.area;
			momentum.comp[2]=flux[3];//-(Pref*face.normal.comp[2])*face.area;
			bc.region[face.bc].momentum-=momentum;
		}

		// Fill in residual (rhs vector)
		for (int i=0;i<nSolVar;++i) {
			row=grid.cell[parent].globalId*nSolVar+i;
			value=flux[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==-1) { // TODO what if a ghost face??
				row=grid.cell[neighbor].globalId*nSolVar+i;
				value*=-1.;
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}

		// TODO flux for the diffusive flux jacobian is slightly different
		// due to perturbation, account for that here
		// Think of a better (cheaper) way to do this
		for (int m=0;m<7;++m) {
			flux[m]=0.;
			face_state_adjust(left,right,face,f,m);
		}
		convective_face_flux(left,right,face,f,flux);
		diffusive_face_flux(left,right,face,f,flux);

		for (int i=0;i<nSolVar;++i) {

			// Flush face fluxes
			for (int m=0;m<7;++m) fluxPlus[m]=0.;
			// Perturb left variable
			left_state_perturb(left,right,face,f,i,epsilon);
			// Adjust face averages and gradients (crude)
			face_state_adjust(left,right,face,f,i);

			convective_face_flux(left,right,face,f,fluxPlus);
			diffusive_face_flux(left,right,face,f,fluxPlus);

			// Restore
			left_state_perturb(left,right,face,f,i,-1.*epsilon);
			face_state_adjust(left,right,face,f,i);

			// Add change of flux (flux Jacobian) to implicit operator
			for (int j=0;j<nSolVar;++j) {
				row=grid.cell[parent].globalId*nSolVar+j;
				col=grid.cell[parent].globalId*nSolVar+i;
				value=(flux[j]-fluxPlus[j])/epsilon;
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				if (face.bc==-1) { // TODO what if a ghost face??
					row=grid.cell[neighbor].globalId*nSolVar+j;
					value*=-1.;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				}
			} // for j

			if (face.bc==-1) {

				// Flush face fluxes
				for (int m=0;m<7;++m) fluxPlus[m]=0.;

				// Perturb right variable
				right_state_perturb(right,face,i,epsilon);
				// Adjust face averages and gradients (crude)
				face_state_adjust(left,right,face,f,i);

				convective_face_flux(left,right,face,f,fluxPlus);
				diffusive_face_flux(left,right,face,f,fluxPlus);

				// Restore
				right_state_perturb(right,face,i,-1.*epsilon);
				face_state_adjust(left,right,face,f,i);

				// Add change of flux (flux Jacobian) to implicit operator
				for (int j=0;j<nSolVar;++j) {
					row=grid.cell[neighbor].globalId*nSolVar+j;
					col=grid.cell[neighbor].globalId*nSolVar+i;
					value=(fluxPlus[j]-flux[j])/epsilon;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					row=grid.cell[parent].globalId*nSolVar+j;
					value*=-1.;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				} // for j
			} // if grid.face[f].bc==-1
		} // for i

	} // for faces

	MatAssemblyBegin(impOP,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FLUSH_ASSEMBLY);

	return;
} // end function


void left_state_update(Cell_State &left,Face_State &face,unsigned int f) {

	unsigned int parent;
	Vec3D deltaV;
	double delta[7];

	parent=grid.face[f].parent;

	if (order=="second") {
		for (unsigned int i=0;i<7;++i) {
			delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].limited_grad[i]);
		}
	} else {
		for (unsigned int i=0;i<7;++i) delta[i]=0.;
	}

	deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
	
	// Set left primitive variables
	left.rho=grid.cell[parent].rho+delta[0];
	left.v=grid.cell[parent].v+deltaV;
	left.p=grid.cell[parent].p+delta[4];
	left.k=grid.cell[parent].k+delta[5];
	left.omega=grid.cell[parent].omega+delta[6];
	left.a=sqrt(Gamma*(left.p+Pref)/left.rho);
	left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);

	left.vN.comp[0]=left.v.dot(face.normal);
	left.vN.comp[1]=left.v.dot(face.tangent1);
	left.vN.comp[2]=left.v.dot(face.tangent2);
	
	return;
} // end left_state_update

void right_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f) {

	if (face.bc==-1) {// internal face

		Vec3D deltaV;
		double delta[7];
		unsigned int neighbor;

  		neighbor=grid.face[f].neighbor;

		if (order=="second") {
			for (unsigned int i=0;i<7;++i) {
				delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].limited_grad[i]);
			}
		} else {
			for (unsigned int i=0;i<7;++i) delta[i]=0.;
		}

		deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];

		// Set right primitive variables
		right.rho=grid.cell[neighbor].rho+delta[0];
		right.v=grid.cell[neighbor].v+deltaV;
		right.p=grid.cell[neighbor].p+delta[4];
		right.k=grid.cell[neighbor].k+delta[5];
		right.omega=grid.cell[neighbor].omega+delta[6];
		
	} else if (face.bc>=0) { // boundary face

		if (bc.region[face.bc].type=="outlet") {
			right.rho=left.rho;
			right.v=left.v;
			right.p=left.p;
			right.k=left.k;
			right.omega=left.omega;
			if (bc.region[face.bc].kind=="fixedPressure") {
				double Mach=(left.v.dot(face.normal))/left.a;
				if (Mach<1.) right.p=bc.region[face.bc].p;
			}
		} else if (bc.region[face.bc].type=="slip" | bc.region[face.bc].type=="symmetry") {
			right.rho=left.rho;
			right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
			right.p=left.p;
			right.k=left.k;
			right.omega=left.omega;
		} else if (bc.region[face.bc].type=="noslip") {
			right.rho=left.rho;
			right.v=-1.*left.v;
			right.p=left.p;
			right.k=left.k;
			right.omega=left.omega;
		} else if (bc.region[face.bc].type=="inlet") {
			double Mach=(left.v.dot(face.normal))/left.a;
			right.rho=bc.region[face.bc].rho;
			right.v=bc.region[face.bc].v;
			right.p=left.p;
			if (Mach<=-1.) right.p=bc.region[face.bc].p;
			right.k=bc.region[face.bc].k;
			right.omega=bc.region[face.bc].omega;
		}

	} else { // partition boundary
		int g=-1*face.bc-3;
		Vec3D deltaV;
		double delta[7];
		if (order=="second") {
			for (unsigned int i=0;i<7;++i) {
				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].limited_grad[i]);
			}
		}

		deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
		right.rho=grid.ghost[g].rho+delta[0];	
		right.v=grid.ghost[g].v+deltaV;
		right.p=grid.ghost[g].p+delta[4];
		right.k=grid.ghost[g].k+delta[5];
		right.omega=grid.ghost[g].omega+delta[6];
	}

	right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
	right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
	right.vN.comp[0]=right.v.dot(face.normal);
	right.vN.comp[1]=right.v.dot(face.tangent1);
	right.vN.comp[2]=right.v.dot(face.tangent2);
	
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
	if (face.bc==-1) { rightCentroid=grid.cell[grid.face[f].neighbor].centroid;}
	else if (face.bc>=0) {
		rightCentroid=leftCentroid+2.*(grid.face[f].centroid-leftCentroid).dot(face.normal)*face.normal;
	}
	else { int g=-1*face.bc-3; rightCentroid=grid.ghost[g].centroid;}
	face.left2right=rightCentroid-leftCentroid;
	return;
} // end face_geom_update


void face_state_update(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f) {

	map<int,double>::iterator fit;
	// Find face averaged variables
	face.gradU=0.; face.gradV=0.; face.gradW=0.; face.gradK=0.; face.gradOmega=0.;
	for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			face.gradU+=(*fit).second*grid.cell[(*fit).first].grad[1];
			face.gradV+=(*fit).second*grid.cell[(*fit).first].grad[2];
			face.gradW+=(*fit).second*grid.cell[(*fit).first].grad[3];
			face.gradK+=(*fit).second*grid.cell[(*fit).first].grad[5];
			face.gradOmega+=(*fit).second*grid.cell[(*fit).first].grad[6];
		} else { // if contribution is coming from a ghost cell
			face.gradU+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[1];
			face.gradV+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[2];
			face.gradW+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[3];
			face.gradK+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[5];
			face.gradOmega+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[6];
		}
	}

	// Boundary conditions are already taken care of in right state update
	face.rho=0.5*(left.rho+right.rho);
	face.v=0.5*(left.v+right.v);
	face.p=0.5*(left.p+right.p);
	face.k=0.5*(left.k+right.k);
	face.omega=0.5*(left.omega+right.omega);
	
	return;
} // end face_state_update

void left_state_perturb(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,int var,double epsilon) {

	if (var==0) { // perturb rho
		left.rho+=epsilon;
		left.a=sqrt(Gamma*(left.p+Pref)/left.rho);
		left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
		if (face.bc>=0 && bc.region[face.bc].type!="inlet") {
			right.rho+=epsilon;
			right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
			right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
		}
	} else if (var==1) { // perturb u
		left.v.comp[0]+=epsilon;
		left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
		left.vN.comp[0]=left.v.dot(face.normal);
		left.vN.comp[1]=left.v.dot(face.tangent1);
		left.vN.comp[2]=left.v.dot(face.tangent2);
		if (face.bc>=0) {
			if (bc.region[face.bc].type=="outlet")  { right.v.comp[0]+=epsilon; }
			else if (bc.region[face.bc].type=="slip" | bc.region[face.bc].type=="symmetry") { right.v=left.v-2.*left.v.dot(face.normal)*face.normal; }
			else if (bc.region[face.bc].type=="noslip") { right.v.comp[0]-=epsilon; }
			right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
			right.vN.comp[0]=right.v.dot(face.normal);
			right.vN.comp[1]=right.v.dot(face.tangent1);
			right.vN.comp[2]=right.v.dot(face.tangent2);
		}
	} else if (var==2) { // perturb v
		left.v.comp[1]+=epsilon;
		left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
		left.vN.comp[0]=left.v.dot(face.normal);
		left.vN.comp[1]=left.v.dot(face.tangent1);
		left.vN.comp[2]=left.v.dot(face.tangent2);
		if (face.bc>=0) {
			if (bc.region[face.bc].type=="outlet")  { right.v.comp[1]+=epsilon; }
			else if (bc.region[face.bc].type=="slip" | bc.region[face.bc].type=="symmetry") { right.v=left.v-2.*left.v.dot(face.normal)*face.normal; }
			else if (bc.region[face.bc].type=="noslip") { right.v.comp[1]-=epsilon; }
			right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
			right.vN.comp[0]=right.v.dot(face.normal);
			right.vN.comp[1]=right.v.dot(face.tangent1);
			right.vN.comp[2]=right.v.dot(face.tangent2);
		}
	} else if (var==3) { // perturb w
		left.v.comp[2]+=epsilon;
		left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
		left.vN.comp[0]=left.v.dot(face.normal);
		left.vN.comp[1]=left.v.dot(face.tangent1);
		left.vN.comp[2]=left.v.dot(face.tangent2);
		if (face.bc>=0) {
			if (bc.region[face.bc].type=="outlet")  { right.v.comp[2]+=epsilon; }
			else if (bc.region[face.bc].type=="slip"  | bc.region[face.bc].type=="symmetry") { right.v=left.v-2.*left.v.dot(face.normal)*face.normal; }
			else if (bc.region[face.bc].type=="noslip") { right.v.comp[2]-=epsilon; }
			right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
			right.vN.comp[0]=right.v.dot(face.normal);
			right.vN.comp[1]=right.v.dot(face.tangent1);
			right.vN.comp[2]=right.v.dot(face.tangent2);
		}
	} else if (var==4) { // perturb p
		left.p+=epsilon;
		left.a=sqrt(Gamma*(left.p+Pref)/left.rho);
		left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
		if (face.bc>=0) {
			if (bc.region[face.bc].type=="outlet") {
				if (bc.region[face.bc].kind=="fixedPressure") {
					double Mach=(left.v.dot(face.normal))/left.a;
					if (Mach<1.) right.p=bc.region[face.bc].p;
				}
				right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
			} else if (bc.region[face.bc].type=="inlet") {
				double Mach=(left.v.dot(face.normal))/left.a;
				right.p=left.p;
				if (Mach<=-1.) right.p=bc.region[face.bc].p;
			} else {
				right.p+=epsilon;
			}
			right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
			right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
		}
	} else if (var==5) { // perturb k
		left.k+=epsilon;
		if (face.bc>=0 && bc.region[face.bc].type!="inlet") right.k+=epsilon;
	} else if (var==6) { // perturb omega
		left.omega+=epsilon;
		if (face.bc>=0 && bc.region[face.bc].type!="inlet") right.omega+=epsilon;
	}
	
	return;
} // end left_state_perturb

void right_state_perturb(Cell_State &right,Face_State &face,int var,double epsilon) {

	// TODO use a case switch here instead of if
	// TODO same goes for left_state_perturb
	if (var==0) { // perturb rho
		right.rho+=epsilon;
		right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
		right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
	} else if (var==1) { // perturb u
		right.v.comp[0]+=epsilon;
		right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
		right.vN.comp[0]=right.v.dot(face.normal);
		right.vN.comp[1]=right.v.dot(face.tangent1);
		right.vN.comp[2]=right.v.dot(face.tangent2);
	} else if (var==2) { // perturb v
		right.v.comp[1]+=epsilon;
		right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
		right.vN.comp[0]=right.v.dot(face.normal);
		right.vN.comp[1]=right.v.dot(face.tangent1);
		right.vN.comp[2]=right.v.dot(face.tangent2);
	} else if (var==3) { // perturb w
		right.v.comp[2]+=epsilon;
		right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
		right.vN.comp[0]=right.v.dot(face.normal);
		right.vN.comp[1]=right.v.dot(face.tangent1);
		right.vN.comp[2]=right.v.dot(face.tangent2);
	} else if (var==4) { // perturb p
		right.p+=epsilon;
		right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
		right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
	} else if (var==5) { // perturb k
		right.k+=epsilon;
	} else if (var==6) { // perturb omega
		right.omega+=epsilon;
	}
	
	return;
} // end right_state_perturb

void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,int var) {

	double distance;
	Vec3D direction;
	
	if (var!=0 && var!=4) {
		distance=fabs(face.left2right);
		direction=face.left2right/distance;
	}

	switch (var)
	{
		case 0 : // rho
			face.rho=0.5*(left.rho+right.rho);
			break;
		case 1 : // u
			face.v.comp[0]=0.5*(left.v.comp[0]+right.v.comp[0]);
			face.gradU-=face.gradU.dot(direction)*direction;
			face.gradU+=((right.v.comp[0]-left.v.comp[0])/distance)*direction;
			break;
		case 2 : // v
			face.v.comp[1]=0.5*(left.v.comp[1]+right.v.comp[1]);
			face.gradV-=face.gradV.dot(direction)*direction;
			face.gradV+=((right.v.comp[1]-left.v.comp[1])/distance)*direction;
			break;
		case 3 : // w
			face.v.comp[2]=0.5*(left.v.comp[2]+right.v.comp[2]);
			face.gradW-=face.gradW.dot(direction)*direction;
			face.gradW+=((right.v.comp[2]-left.v.comp[2])/distance)*direction;
			break;
		case 4 : // p
			face.p=0.5*(left.p+right.p);
			break;
		case 5 : // k
			face.k=0.5*(left.k+right.k);
			face.gradK-=face.gradK.dot(direction)*direction;
			face.gradK+=((right.k-left.k)/distance)*direction;
			break;
		case 6 : // omega
			face.omega=0.5*(left.omega+right.omega);
			face.gradOmega-=face.gradOmega.dot(direction)*direction;
			face.gradOmega+=((right.omega-left.omega)/distance)*direction;
			break;
	}

	return;
} // end face_state_adjust
