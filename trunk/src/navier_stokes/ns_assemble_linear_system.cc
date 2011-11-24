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

double order_factor;

namespace ns_state {
	NS_Cell_State left,right,leftPlus,rightPlus;
	NS_Face_State face;
	NS_Fluxes flux, fluxPlus;
	vector<double> jacobianLeft,jacobianRight,sourceJacLeft,sourceJacRight;
	vector<double> sourceLeft,sourceRight,sourceLeftPlus,sourceRightPlus;
	bool doLeftSourceJac,doRightSourceJac;
}

void NavierStokes::assemble_linear_system(void) {

	MatZeroEntries(impOP);
	
	using namespace ns_state;
	using ns_state::left;
	using ns_state::right;
		
	int parent,neighbor,f;
	int row,col;
	
	vector<bool> cellVisited;
	for (int c=0;c<grid[gid].cellCount;++c) cellVisited.push_back(false);
	
	small_number=10.*sqrt(std::numeric_limits<double>::epsilon());
	
	PetscScalar value;

	flux.convective.resize(5);
	flux.diffusive.resize(5);
	fluxPlus.convective.resize(5);
	fluxPlus.diffusive.resize(5);
	jacobianLeft.resize(5);
	jacobianRight.resize(5);
	left.update.resize(5);
	right.update.resize(5);
	leftPlus.update.resize(5);
	rightPlus.update.resize(5);
	sourceLeft.resize(5);
	sourceRight.resize(5);
	sourceLeftPlus.resize(5);
	sourceRightPlus.resize(5);
	sourceJacLeft.resize(5);
	sourceJacRight.resize(5);
	
	for (int m=0;m<5;++m) flux.diffusive[m]=0.;
	
	bool implicit=true;
	// if (TIME_INTEGRATOR==FORWARD_EULER) implicit=false; // TODO : properly implement explicit solver

	// Loop through faces
	for (f=0;f<grid[gid].faceCount;++f) {
		
		order_factor=1.; // This is used to kill gradients during Jacobian calculation if first order option is selected
		
		doLeftSourceJac=false; doRightSourceJac=false;
		for (int m=0;m<5;++m) { 
			sourceLeft[m]=0.;
			sourceRight[m]=0.;
		}
		parent=grid[gid].face[f].parent; neighbor=grid[gid].face[f].neighbor;

		// Populate the state caches
		face_geom_update(face,f);
		left_state_update(left,face);
		right_state_update(left,right,face);
		face_state_update(left,right,face);
		// Get unperturbed flux values
		convective_face_flux(left,right,face,&flux.convective[0]);
		diffusive_face_flux(left,right,face,&flux.diffusive[0]);
		
		// Add Sources
		if (!cellVisited[parent]){
			sources(left,&sourceLeft[0]);
			cellVisited[parent]=true;
			doLeftSourceJac=true;
		}

		if(face.bc==INTERNAL_FACE && !cellVisited[neighbor]) {
			sources(right,&sourceRight[0]);
			cellVisited[neighbor]=true;
			doRightSourceJac=true;
		}
	
		// Integrate boundary loads
		if ((timeStep) % loads[gid].frequency == 0) {
			Vec3D temp;
			for (int b=0;b<loads[gid].include_bcs.size();++b) {
				if (face.bc==loads[gid].include_bcs[b]) {
					for (int i=0;i<3;++i) temp[i]=flux.convective[i+1]-flux.diffusive[i+1];
					loads[gid].force[b]+=temp;
					loads[gid].moment[b]+=(grid[gid].face[face.index].centroid-loads[gid].moment_center).cross(temp);
					break;
				}
			}
		}
	
		// Fill in surface information
		
		if (face.bc>=0) {
			if (!qdot.fixedonBC[face.bc]) qdot.bc(face.bc,face.index)=(flux.convective[4]-flux.diffusive[4])/face.area;
			if (!mdot.fixedonBC[face.bc]) mdot.bc(face.bc,face.index)=(flux.convective[0]-flux.diffusive[0])/face.area;
			for (int i=0;i<3;++i) tau.bc(face.bc,face.index)[i]=-flux.diffusive[i+1]/face.area;
		}
		
		// Fill in rhs vector
		for (int i=0;i<5;++i) {
			row=(grid[gid].myOffset+parent)*5+i;
			value=flux.diffusive[i]-flux.convective[i]+sourceLeft[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid[gid].face[f].bc==INTERNAL_FACE) { 
				row=(grid[gid].myOffset+neighbor)*5+i;
				value=-1.*(flux.diffusive[i]-flux.convective[i])+sourceRight[i];
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}

		//if (implicit && ps_timeStep==1) { // TODO: Get this working

			for (int i=0;i<5;++i) { // perturb each variable
				
				for (int m=0;m<5;++m) { 
					sourceJacLeft[m]=0.;
					sourceJacRight[m]=0.;
				}
				
				get_jacobians(i);
	
				// Add change of flux (flux Jacobian) to implicit operator
				for (int j=0;j<5;++j) {
					col=(grid[gid].myOffset+parent)*5+i; // Effect of parent ith var perturbation
					row=(grid[gid].myOffset+parent)*5+j; // on parent jth flux
					value=-1.*jacobianLeft[j];
					//if (doLeftSourceJac) value-=sourceJacLeft[j];
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					if (face.bc==INTERNAL_FACE) { 
						row=(grid[gid].myOffset+neighbor)*5+j; // on neighbor jth flux
						value=jacobianLeft[j]; 
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					}
				} // for j
	
				if (face.bc==INTERNAL_FACE || face.bc==GHOST_FACE) { 
	
					// Add change of flux (flux Jacobian) to implicit operator
					for (int j=0;j<5;++j) {
						if (face.bc==INTERNAL_FACE) { 
							col=(grid[gid].myOffset+neighbor)*5+i; // Effect of neighbor ith var perturbation
							row=(grid[gid].myOffset+neighbor)*5+j; // on neighbor jth flux
							value=jacobianRight[j];
							//if (doRightSourceJac) value-=sourceJacRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
							row=(grid[gid].myOffset+parent)*5+j; // on parent jth flux
							value=-1.*jacobianRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						} else  { // Ghost (only add effect on parent cell, effect on itself is taken care of in its own partition
							row=(grid[gid].myOffset+parent)*5+j;
							col=(grid[gid].ghost[-1*neighbor-1].matrix_id)*5+i;
							value=-1.*jacobianRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						}
					}
						
				} // if 
				
			} // for i (each perturbed variable)
		//} // if implicit

	} // for faces
	
	return;
} // end function

void NavierStokes::get_jacobians(const int var) {

	using namespace ns_state;
	using ns_state::left;
	using ns_state::right;
	double epsilon;
	double factor=0.001;
	
	if (jac_order==FIRST) order_factor=0.;
	
	leftPlus=left;
	rightPlus=right;
	
	for (int m=0;m<5;++m) {
		sourceLeftPlus[m]=sourceLeft[m];
		sourceRightPlus[m]=sourceRight[m];
		sourceJacLeft[m]=0.;
		sourceJacRight[m]=0.;
	}
	
	if (left.update[var]>small_number) {epsilon=max(small_number,factor*left.update[var]);}
	else if (left.update[var]<-small_number) {epsilon=min(-small_number,factor*left.update[var]); }
	else {epsilon=small_number;}
	
	// Perturb left state
	state_perturb(leftPlus,face,var,epsilon);
	// If right state is a boundary, correct the condition according to changes in left state
	if (face.bc>=0) right_state_update(leftPlus,rightPlus,face);

	convective_face_flux(leftPlus,rightPlus,face,&fluxPlus.convective[0]);
	
	face_state_adjust(leftPlus,rightPlus,face,var);
	diffusive_face_flux(leftPlus,rightPlus,face,&fluxPlus.diffusive[0]);
	
	if (doLeftSourceJac) {
		sources(left,&sourceLeft[0],true);
		sources(leftPlus,&sourceLeftPlus[0],true);
	}
	
	for (int j=0;j<5;++j){
		jacobianLeft[j]=(fluxPlus.diffusive[j]-flux.diffusive[j]-fluxPlus.convective[j]+flux.convective[j])/epsilon;
		if (doLeftSourceJac) sourceJacLeft[j]=(sourceLeftPlus[j]-sourceLeft[j])/epsilon;
	}
	
	if (face.bc==INTERNAL_FACE || face.bc==GHOST_FACE) { 
		
		if (right.update[var]>small_number) {epsilon=max(small_number,factor*right.update[var]);}
		else if (right.update[var]<-small_number) {epsilon=min(-small_number,factor*right.update[var]); }
		else {epsilon=small_number;}
				
		state_perturb(rightPlus,face,var,epsilon);
		
		convective_face_flux(left,rightPlus,face,&fluxPlus.convective[0]);
		face_state_adjust(left,rightPlus,face,var);
		diffusive_face_flux(left,rightPlus,face,&fluxPlus.diffusive[0]);
		
		if (doRightSourceJac) {
			sources(right,&sourceRight[0],true);
			sources(rightPlus,&sourceRightPlus[0],true);
		}
		
		for (int j=0;j<5;++j){
			jacobianRight[j]=(fluxPlus.diffusive[j]-flux.diffusive[j]-fluxPlus.convective[j]+flux.convective[j])/epsilon;
			if (doRightSourceJac) sourceJacRight[j]=(sourceRightPlus[j]-sourceRight[j])/epsilon;
		}
	
	}
	
	face_state_adjust(left,right,face,var);
	
}

void NavierStokes::left_state_update(NS_Cell_State &left,NS_Face_State &face) {
	
	int parent=grid[gid].face[face.index].parent;

	Vec3D cell2face=order_factor*grid[gid].face[face.index].centroid-grid[gid].cell[parent].centroid;
	Vec3D deltaV;
	left.p_center=p.cell(parent);
	left.p=left.p_center+limiter[0].cell(parent)*cell2face.dot(gradp.cell(parent));
	deltaV[0]=limiter[1].cell(parent)*cell2face.dot(gradu.cell(parent));
	deltaV[1]=limiter[2].cell(parent)*cell2face.dot(gradv.cell(parent));
	deltaV[2]=limiter[3].cell(parent)*cell2face.dot(gradw.cell(parent));
	left.V_center=V.cell(parent);
	left.V=left.V_center+deltaV;
	left.T_center=T.cell(parent);
	left.T=left.T_center+limiter[4].cell(parent)*cell2face.dot(gradT.cell(parent));
	left.rho=material.rho(left.p,left.T);
	
	for (int i=0;i<5;++i) left.update[i]=update[i].cell(parent);
		
	left.a=material.a(left.p,left.T);
	left.H=left.a*left.a/(material.gamma-1.)+0.5*left.V.dot(left.V);
	left.Vn[0]=left.V.dot(face.normal);
	left.Vn[1]=left.V.dot(face.tangent1);
	left.Vn[2]=left.V.dot(face.tangent2);
	left.volume=grid[gid].cell[parent].volume;
	
	return;
}

void NavierStokes::right_state_update(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {

	if (face.bc==INTERNAL_FACE) {// internal face

		int neighbor=grid[gid].face[face.index].neighbor;
		Vec3D cell2face=order_factor*grid[gid].face[face.index].centroid-grid[gid].cell[neighbor].centroid;
		Vec3D deltaV;
		right.p_center=p.cell(neighbor);
		right.p=right.p_center+limiter[0].cell(neighbor)*cell2face.dot(gradp.cell(neighbor));
		deltaV[0]=limiter[1].cell(neighbor)*cell2face.dot(gradu.cell(neighbor));
		deltaV[1]=limiter[2].cell(neighbor)*cell2face.dot(gradv.cell(neighbor));
		deltaV[2]=limiter[3].cell(neighbor)*cell2face.dot(gradw.cell(neighbor));
		right.V_center=V.cell(neighbor);
		right.V=right.V_center+deltaV;
		right.T_center=T.cell(neighbor);
		right.T=right.T_center+limiter[4].cell(neighbor)*cell2face.dot(gradT.cell(neighbor));
		right.rho=material.rho(right.p,right.T);
		right.volume=grid[gid].cell[neighbor].volume;
		
		for (int i=0;i<5;++i) right.update[i]=update[i].cell(neighbor);
		
	} else if (face.bc>=0) { // boundary face

		for (int i=0;i<5;++i) right.update[i]=0.;
		apply_bcs(left,right,face);
		
	} else { // partition boundary

		int g=-1*grid[gid].face[face.index].neighbor-1; // ghost cell index
		
		Vec3D cell2face=order_factor*grid[gid].face[face.index].centroid-grid[gid].ghost[g].centroid;
		Vec3D deltaV;
		right.p_center=p.ghost(g);
		right.p=right.p_center+limiter[0].ghost(g)*cell2face.dot(gradp.ghost(g));
		deltaV[0]=limiter[1].ghost(g)*cell2face.dot(gradu.ghost(g));
		deltaV[1]=limiter[2].ghost(g)*cell2face.dot(gradv.ghost(g));
		deltaV[2]=limiter[3].ghost(g)*cell2face.dot(gradw.ghost(g));
		right.V_center=V.ghost(g);
		right.V=right.V_center+deltaV;
		right.T_center=T.ghost(g);
		right.T=right.T_center+limiter[4].ghost(g)*cell2face.dot(gradT.ghost(g));
		right.rho=material.rho(right.p,right.T);
		right.volume=grid[gid].ghost[g].volume;
		
		for (int i=0;i<5;++i) right.update[i]=update[i].ghost(g);
		
	}
	
	right.a=material.a(right.p,right.T);
	right.H=right.a*right.a/(material.gamma-1.)+0.5*right.V.dot(right.V);
	right.Vn[0]=right.V.dot(face.normal);
	right.Vn[1]=right.V.dot(face.tangent1);
	right.Vn[2]=right.V.dot(face.tangent2);
		
	return;
}

void NavierStokes::face_geom_update(NS_Face_State &face,int f) {
	face.index=f;
	face.normal=grid[gid].face[f].normal;
	if (grid[gid].face[f].nodeCount==4) {
		face.tangent1=((grid[gid].faceNode(f,0)+grid[gid].faceNode(f,1))-(grid[gid].faceNode(f,2)+grid[gid].faceNode(f,3))).norm();
	} else {
		face.tangent1=(0.5*(grid[gid].faceNode(f,0)+grid[gid].faceNode(f,1))-grid[gid].face[f].centroid).norm();
	}
	// Cross the tangent vector with the normal vector to get the second tangent
	face.tangent2=((face.normal).cross(face.tangent1)).norm();
	face.area=grid[gid].face[f].area;
	face.bc=grid[gid].face[f].bc;
	Vec3D leftCentroid,rightCentroid;
	leftCentroid=grid[gid].cell[grid[gid].face[f].parent].centroid;
	if (face.bc==INTERNAL_FACE) { rightCentroid=grid[gid].cell[grid[gid].face[f].neighbor].centroid;}
	else if (face.bc>=0) {
		if (bc[gid][face.bc].type==OUTLET || bc[gid][face.bc].type==INLET) {
			rightCentroid=leftCentroid+2.*(grid[gid].face[f].centroid-leftCentroid);
		} else {
			rightCentroid=leftCentroid+2.*(grid[gid].face[f].centroid-leftCentroid).dot(face.normal)*face.normal;
		}
	} else { rightCentroid=grid[gid].ghost[-1*grid[gid].face[f].neighbor-1].centroid;}
	face.left2right=rightCentroid-leftCentroid;
	return;
} // end face_geom_update

void NavierStokes::face_state_update(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {

	face.gradu=gradu.face(face.index);
	face.gradv=gradv.face(face.index);
	face.gradw=gradw.face(face.index);
	face.gradT=gradT.face(face.index);
	
	//  The old way of doing deferred correction (REMOVE after enough testing)
	//face.gradu-=face.gradu.dot(face.normal)*face.normal;
	//face.gradu+=((right.V_center[0]-left.V_center[0])/(face.left2right.dot(face.normal)))*face.normal;

	Vec3D l2rnormal=face.left2right;
	l2rnormal=l2rnormal.norm();
	double l2rmag=fabs(face.left2right);

	face.gradu-=face.gradu.dot(l2rnormal)*l2rnormal;
	face.gradu+=((right.V_center[0]-left.V_center[0])/(l2rmag))*l2rnormal;
	
	face.gradv-=face.gradv.dot(l2rnormal)*l2rnormal;
	face.gradv+=((right.V_center[1]-left.V_center[1])/(l2rmag))*l2rnormal;
	
	face.gradw-=face.gradw.dot(l2rnormal)*l2rnormal;
	face.gradw+=((right.V_center[2]-left.V_center[2])/(l2rmag))*l2rnormal;
	
	face.gradT-=face.gradT.dot(l2rnormal)*l2rnormal;
	face.gradT+=((right.T_center-left.T_center)/(l2rmag))*l2rnormal;
	 
	// Boundary conditions are already taken care of in right state update

	// TODO: In first order, this won't be a good averaging
	
	//if (face.bc>=0) {
		face.p=0.5*(left.p+right.p);
		face.V=0.5*(left.V+right.V);
		face.T=0.5*(left.T+right.T);
	//} else {
	//	face.p=p.face(face.index);
	//	face.V=V.face(face.index);
	//	face.T=T.face(face.index);
	//}
	face.mu=material.viscosity(face.T);
	face.lambda=material.therm_cond(face.T);
	
	return;
} // end face_state_update

void NavierStokes::state_perturb(NS_Cell_State &state,NS_Face_State &face,int var,double epsilon) {
	
	switch (var)
	{
		case 0 : // p
			state.p+=epsilon;
			state.rho=material.rho(state.p,state.T);
			state.a=material.a(state.p,state.T);
			break;
		case 1 : // u
			state.V[0]+=epsilon;
			state.V_center[0]+=epsilon;
			state.Vn[0]=state.V.dot(face.normal);
			state.Vn[1]=state.V.dot(face.tangent1);
			state.Vn[2]=state.V.dot(face.tangent2);
			break;
		case 2 : // v
			state.V[1]+=epsilon;
			state.V_center[1]+=epsilon;
			state.Vn[0]=state.V.dot(face.normal);
			state.Vn[1]=state.V.dot(face.tangent1);
			state.Vn[2]=state.V.dot(face.tangent2);
			break;
		case 3 : // w
			state.V[2]+=epsilon;
			state.V_center[2]+=epsilon;
			state.Vn[0]=state.V.dot(face.normal);
			state.Vn[1]=state.V.dot(face.tangent1);
			state.Vn[2]=state.V.dot(face.tangent2);
			break;
		case 4 : // T
			state.T+=epsilon;
			state.T_center+=epsilon;
			state.rho=material.rho(state.p,state.T);
			state.a=material.a(state.p,state.T);
			break;
 	}
	
	state.H=state.a*state.a/(material.gamma-1.)+0.5*state.V.dot(state.V);
	
	return;
} // end state_perturb

void NavierStokes::face_state_adjust(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,int var) {


	Vec3D l2rnormal=face.left2right;
	l2rnormal=l2rnormal.norm();
	double l2rmag=fabs(face.left2right);

	switch (var)
	{
		case 0 : // p
			face.p=0.5*(left.p+right.p);
			//if (face.bc<0) face.p=p.face(face.index);
			break;
		case 1 : // u
			face.V[0]=0.5*(left.V[0]+right.V[0]);
			//if (face.bc<0) face.V[0]=V.face(face.index)[0];
			face.gradu-=face.gradu.dot(l2rnormal)*l2rnormal;
			face.gradu+=((right.V_center[0]-left.V_center[0])/(l2rmag))*l2rnormal;
		case 2 : // v
			face.V[1]=0.5*(left.V[1]+right.V[1]);
			//if (face.bc<0) face.V[1]=V.face(face.index)[1];
			face.gradv-=face.gradv.dot(l2rnormal)*l2rnormal;
			face.gradv+=((right.V_center[1]-left.V_center[1])/(l2rmag))*l2rnormal;
			break;
		case 3 : // w
			face.V[2]=0.5*(left.V[2]+right.V[2]);
			//if (face.bc<0) face.V[2]=V.face(face.index)[2];
			face.gradw-=face.gradw.dot(l2rnormal)*l2rnormal;
			face.gradw+=((right.V_center[2]-left.V_center[2])/(l2rmag))*l2rnormal;
			break;
		case 4 : // T
			face.T=0.5*(left.T+right.T);
			//if (face.bc<0) face.T=T.face(face.index);
			face.gradT-=face.gradT.dot(l2rnormal)*l2rnormal;
			face.gradT+=((right.T_center-left.T_center)/(l2rmag))*l2rnormal;
			face.mu=material.viscosity(face.T);
			face.lambda=material.therm_cond(face.T);
			break;
	}

	return;
} // end face_state_adjust
