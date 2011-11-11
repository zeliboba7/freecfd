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
#include "rans.h"

extern vector<RANS> rans;

NavierStokes::NavierStokes (void) {
	// Empty constructor
	return;
}

void NavierStokes::initialize (int ps_max) {
	
	ps_step_max=ps_max;
	nVars=5;
	rtol=input.section("grid",gid).subsection("navierstokes").get_double("relativetolerance");
	abstol=input.section("grid",gid).subsection("navierstokes").get_double("absolutetolerance");
	maxits=input.section("grid",gid).subsection("navierstokes").get_int("maximumiterations");

	if (input.section("grid",gid).subsection("navierstokes").get_string("limiter")=="vk") {
		limiter_function=VK;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("limiter")=="bj") {
		limiter_function=BJ;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("limiter")=="minmod") {
		limiter_function=MINMOD;
	} else {
		limiter_function=NONE;
	}
	limiter_threshold=input.section("grid",gid).subsection("navierstokes").get_double("limiterthreshold");
	
	if (input.section("grid",gid).subsection("navierstokes").get_string("order")=="first") {
		order=FIRST;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("order")=="second") {
		order=SECOND;
	}
	
	if (input.section("grid",gid).subsection("navierstokes").get_string("jacobianorder")=="first") {
		jac_order=FIRST;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("jacobianorder")=="second") {
		jac_order=order;
	}
	
	if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="AUSM+up") {
		convective_flux_function=AUSM_PLUS_UP;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="Roe") {
		convective_flux_function=ROE;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="vanLeer") {
		convective_flux_function=VAN_LEER;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="SD-SLAU") {
		convective_flux_function=SD_SLAU;
	}
	
	Minf=input.section("reference").get_double("Mach");
	
	if (input.section("pseudotime").get_string("preconditioner")=="none") {
		preconditioner=NONE;
	} else if (input.section("pseudotime").get_string("preconditioner")=="ws95") {
		preconditioner=WS95;
	}
	
	mpi_init();
	material.set(gid);
	create_vars();
	apply_initial_conditions();
	set_bcs();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	calc_limiter();
	petsc_init();
	first_residuals.resize(3);
	first_ps_residuals.resize(3);
	return;
}

void NavierStokes::solve (int ts,int pts) {
	timeStep=ts;
	ps_step=pts;
	assemble_linear_system();
	time_terms();
	petsc_solve();
	if (turbulent[gid]) rans[gid].solve(timeStep,ps_step);
	update_variables();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	calc_limiter();

	return;
}

void NavierStokes::create_vars (void) {
	// Allocate variables
	// Default option is to store on cell centers and ghosts only
	// Can override by for example: rho.nodeStore=true
	rho.allocate(gid);
	p.allocate(gid);
	T.allocate(gid);
	V.allocate(gid);
	gradrho.allocate(gid);
	gradp.allocate(gid);
	gradT.allocate(gid);
	gradu.allocate(gid);
	gradv.allocate(gid);
	gradw.allocate(gid);
	
	update.resize(5);
	limiter.resize(5);
	for (int i=0; i<5; ++i) {
		update[i].allocate(gid);
		limiter[i].allocate(gid);
	}

	qdot.cellStore=false; qdot.allocate(gid); // is not stored anywhere but the BC
	tau.cellStore=false; tau.allocate(gid);
	mdot.cellStore=false; mdot.faceStore=true; mdot.allocate(gid);	
	
	// This one is for other solvers to find out contribution from either left or right side
	weightL.cellStore=false; weightL.faceStore=true; weightL.allocate(gid);
	
	p_total.cellStore=false; p_total.allocate(gid);	
	T_total.cellStore=false; T_total.allocate(gid);
	
	return;
}

void NavierStokes::apply_initial_conditions (void) {
	// Loop through each initial condition region and apply sequentially
	int count=input.section("grid",gid).subsection("IC",0).count;
	for (int ic=0;ic<count;++ic) {
		// Store the reference to current IC region
		Subsection &region=input.section("grid",gid).subsection("IC",ic);
		
		// Check if gradient test option is chosen
		gradient_test=NONE;
		if (region.get_string("gradienttest")=="linear") {
			gradient_test=LINEAR;
			for (int c=0;c<grid[gid].cellCount;++c) p.cell(c)=grid[gid].cell[c].centroid[0];
		} if (region.get_string("gradienttest")=="quadratic") {
			gradient_test=QUADRATIC;
			// Avoid zero grad points
			// Find min x-coord
			min_x=1.e20;
			max_x=-1.e20;
			for (int c=0;c<grid[gid].cellCount;++c) {
				min_x=min(min_x,grid[gid].cell[c].centroid[0]);
				max_x=max(max_x,grid[gid].cell[c].centroid[0]);
			}
			MPI_Allreduce(&min_x,&min_x,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
			MPI_Allreduce(&max_x,&max_x,1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			
			for (int c=0;c<grid[gid].cellCount;++c) {
				double xc=grid[gid].cell[c].centroid[0];
				p.cell(c)=xc*xc+3.*(max_x-min_x)*xc;
			}
		}
		// If in gradient test mode, break the IC loop
		if (gradient_test!=NONE) break;
		
		
		Vec3D regionV=region.get_Vec3D("V");
		double regionRho,regionP,regionT;
		// Assign specified values
		if (region.get_double("p").is_found) { // If p is specified in the input file
			regionP=region.get_double("p");
			if (region.get_double("rho").is_found) { // Rho is specified in the input file
				//  Check if Temperature is also specified
				if (region.get_double("T").is_found) {
					cerr << "Both p, rho and T can't be specified in initial condition IC_" << ic+1 << " of grid " << gid+1 << endl;
					exit(1);
				}
				regionRho=region.get_double("rho");
				regionT=material.T(regionP,regionRho);
			} else if (region.get_double("T").is_found) {
				regionT=region.get_double("T");
				regionRho=material.rho(regionP,regionT);	
			} else {
				cerr << "Need to specify either rho or T in initial condition IC_" << ic+1 << endl;
				exit(1);
			}
		} else if (region.get_double("rho").is_found) {
			regionRho=region.get_double("rho");
			if (region.get_double("T").is_found) {
				regionT=region.get_double("T");
				regionP=material.p(regionRho,regionT);
			} else {
				cerr << "Need to specify T in initial condition IC_" << ic+1 << endl;
				exit(1);
			}
		} else {
			cerr << "Need to specify either p or rho in initial condition IC_" << ic+1 << endl;
			exit(1);
		}
		
		// If region is specified with a box method
		if (region.get_string("region")=="box") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the box region
				if (withinBox(grid[gid].cell[c].centroid,region.get_Vec3D("corner_1"),region.get_Vec3D("corner_2"))) {
					// Assign specified values
					p.cell(c)=regionP;
					T.cell(c)=regionT;
					rho.cell(c)=regionRho;
					V.cell(c)=regionV;
				}
			}
		} else if (region.get_string("region")=="cylinder") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the cylinder region
				Vec3D axisDirection=region.get_Vec3D("axisdirection");
				axisDirection=axisDirection.norm();
				if (withinCylinder(grid[gid].cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"),axisDirection,region.get_double("height"))) {
					p.cell(c)=regionP;
					T.cell(c)=regionT;
					rho.cell(c)=regionRho;
					// first component of the specified velocity is interpreted as the axial velocity
					// second component of the specified velocity is interpreted as the radial velocity
					// third component of the specified velocity is interpreted as the circumferential velocity
					Vec3D radialDirection=grid[gid].cell[c].centroid-region.get_Vec3D("center");
					radialDirection=(radialDirection-radialDirection.dot(axisDirection)*axisDirection).norm();	
					V.cell(c)=		
					regionV[0]*axisDirection /* axial component */
					+regionV[1]*radialDirection /* radial component */
					+regionV[2]*radialDirection.cross(axisDirection); /* circumferential component */
				}
			}
		} else if (region.get_string("region")=="sphere") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the sphere region
				if (withinSphere(grid[gid].cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"))) {
					p.cell(c)=regionP;
					T.cell(c)=regionT;
					rho.cell(c)=regionRho;
					// first component of the specified velocity is interpreted as the radial velocity
					// second and third components of the specified velocity are ignored
					V.cell(c)=regionV[0]*(grid[gid].cell[c].centroid-region.get_Vec3D("center"));
				}
			}
		}
		
		for (int c=0;c<grid[gid].cellCount;++c) {
			double delta=region.get_double("BLthickness");
			if (delta>0.) {
				double factor=min(1.,grid[gid].cell[c].closest_wall_distance/delta);
				V.cell(c)*=factor;
			}
		}
		
	}
	
	// initialize updates and limiter
	for (int c=0;c<grid[gid].cellCount;++c) {
		for (int i=0;i<5;++i) {
			update[i].cell(c)=0.;
			if (limiter_function==NONE) limiter[i].cell(c)=1.;
			if (order==FIRST) limiter[i].cell(c)=0.;
		}
	}
	for (int g=0;g<grid[gid].ghostCount;++g) {
		for (int i=0;i<5;++i) {
			update[i].ghost(g)=0.;
			if (limiter_function==NONE) limiter[i].ghost(g)=1.;
			if (order==FIRST) limiter[i].ghost(g)=0.;
		}
	}
	return;
}

void NavierStokes::calc_cell_grads (void) {
	vector<Vec3D> grad (3,0.);
	for (int c=0;c<grid[gid].cellCount;++c) {
		gradrho.cell(c)=rho.cell_gradient(c);
		gradp.cell(c)=p.cell_gradient(c);
		gradT.cell(c)=T.cell_gradient(c);

		grad=V.cell_gradient(c);
		gradu.cell(c)[0]=grad[0][0];
		gradu.cell(c)[1]=grad[1][0];
		gradu.cell(c)[2]=grad[2][0];
		gradv.cell(c)[0]=grad[0][1];
		gradv.cell(c)[1]=grad[1][1];
		gradv.cell(c)[2]=grad[2][1];
		gradw.cell(c)[0]=grad[0][2];
		gradw.cell(c)[1]=grad[1][2];
		gradw.cell(c)[2]=grad[2][2];
		
	 }
	return;
}

void NavierStokes::update_variables(void) {
	
	double residuals[3],ps_residuals[3],totalResiduals[3],total_ps_residuals[3];
	for (int i=0;i<3;++i) {
		residuals[i]=0.;
		ps_residuals[i]=0.;
	}

	PetscInt row;
	
	double dt2,dtau2;
	for (int c=0;c<grid[gid].cellCount;++c) {
		for (int i=0;i<5;++i) {
			if (isnan(update[i].cell(c)) || isinf(update[i].cell(c))) {
				cerr << "[E] Divergence detected!...exiting" << endl;
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
		p.cell(c)+=update[0].cell(c);
		T.cell(c)+=update[4].cell(c);
		rho.cell(c)=material.rho(p.cell(c),T.cell(c));
		V.cell(c)[0]+=update[1].cell(c);
		V.cell(c)[1]+=update[2].cell(c);
		V.cell(c)[2]+=update[3].cell(c);
		
		if (ps_step_max>1) {
			dtau2=dtau[gid].cell(c)*dtau[gid].cell(c);
			ps_residuals[0]+=update[0].cell(c)*update[0].cell(c)/dtau2;
			ps_residuals[1]+=update[1].cell(c)*update[1].cell(c)/dtau2+update[2].cell(c)*update[2].cell(c)/dtau2+update[3].cell(c)*update[3].cell(c)/dtau2;
			ps_residuals[2]+=update[4].cell(c)*update[4].cell(c)/dtau2;
		}
		
		// If the last pseudo time step
		if (ps_step_max>1) { // If pseudo time iterations are active
			for (int i=0;i<5;++i) {
				row=(grid[gid].myOffset+c)*5+i;
				VecGetValues(pseudo_delta,1,&row,&update[i].cell(c));
			}
		}
		dt2=dt[gid].cell(c)*dt[gid].cell(c);
		residuals[0]+=update[0].cell(c)*update[0].cell(c)/dt2;
		residuals[1]+=update[1].cell(c)*update[1].cell(c)/dt2+update[2].cell(c)*update[2].cell(c)/dt2+update[3].cell(c)*update[3].cell(c)/dt2;
		residuals[2]+=update[4].cell(c)*update[4].cell(c)/dt2;


		
	} // cell loop
	
	MPI_Allreduce(&residuals,&totalResiduals,3, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if (ps_step_max>1) MPI_Allreduce(&ps_residuals,&total_ps_residuals,3, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		
	if (timeStep==1) for (int i=0;i<3;++i) first_residuals[i]=sqrt(totalResiduals[i]);
	if (ps_step_max>1 && ps_step==1) for (int i=0;i<3;++i) first_ps_residuals[i]=sqrt(total_ps_residuals[i]);
	
	res=0.;
	ps_res=0.;
	for (int i=0;i<3;++i) res+=sqrt(totalResiduals[i])/first_residuals[i]/3.;
	if (ps_step_max>1) for (int i=0;i<3;++i) ps_res+=sqrt(total_ps_residuals[i])/first_ps_residuals[i]/3.;
	
	return;
}

void NavierStokes::find_min_max (void) {

	// Loop the cells and find the local max and min of each primitive vars
	
	qmax[0]=qmin[0]=p.cell(0);
	qmax[1]=qmin[1]=V.cell(0)[0];
	qmax[2]=qmin[2]=V.cell(0)[1];
	qmax[3]=qmin[3]=V.cell(0)[2];
	qmax[4]=qmin[4]=T.cell(0);

	double temp[5];

	for (int c=0;c<grid[gid].cellCount;++c) {
		qmax[0]=max(qmax[0],p.cell(c));		
		qmax[1]=max(qmax[1],V.cell(c)[0]);		
		qmax[2]=max(qmax[2],V.cell(c)[1]);		
		qmax[3]=max(qmax[4],V.cell(c)[2]);		
		qmax[4]=max(qmax[4],p.cell(c));		
		qmin[0]=min(qmin[0],p.cell(c));		
		qmin[1]=min(qmin[1],V.cell(c)[0]);		
		qmin[2]=min(qmin[2],V.cell(c)[1]);		
		qmin[3]=min(qmin[4],V.cell(c)[2]);		
		qmin[4]=min(qmin[4],p.cell(c));		
	}

	MPI_Allreduce(qmax,temp,5, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	for (int i=0;i<5;++i) qmax[i]=temp[i];
	MPI_Allreduce(qmin,temp,5, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	for (int i=0;i<5;++i) qmin[i]=temp[i];

	// DEBUG
	//cout << "MAX\t" << Rank << "\t" << qmax[0] << "\t" << qmax[1] << "\t" << qmax[4] << endl;
	//cout << "MIN\t" << Rank << "\t" << qmin[0] << "\t" << qmin[1] << "\t" << qmin[4] << endl;

	return;
}
