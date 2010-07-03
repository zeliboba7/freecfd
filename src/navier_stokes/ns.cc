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

NavierStokes::NavierStokes (void) {
	// Empty constructor
	return;
}

void NavierStokes::initialize (void) {
	
	rtol=input.section("grid",gid).subsection("navierstokes").get_double("relativetolerance");
	abstol=input.section("grid",gid).subsection("navierstokes").get_double("absolutetolerance");
	maxits=input.section("grid",gid).subsection("navierstokes").get_int("maximumiterations");

	if (input.section("grid",gid).subsection("navierstokes").get_string("limiter")=="venkatakrishnan") {
		limiter_function=VK;
	} else {
		limiter_function=NONE;
	}
	limiter_threshold=input.section("grid",gid).subsection("navierstokes").get_double("limiter_threshold");
	
	if (input.section("grid",gid).subsection("navierstokes").get_string("order")=="first") {
		order=FIRST;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("order")=="second") {
		order=SECOND;
	}
	
	if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="AUSM+up") {
		convective_flux_function=AUSM_PLUS_UP;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="Roe") {
		convective_flux_function=ROE;
	} else if (input.section("grid",gid).subsection("navierstokes").get_string("convectiveflux")=="SD-SLAU") {
		convective_flux_function=SD_SLAU;
	}
	
	Minf=input.section("reference").get_double("Mach");
	
	mpi_init();
	material.set(gid);
	create_vars();
	apply_initial_conditions();
	set_bcs();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	calc_limiter();
	mpi_update_ghost_gradients(); // Needs to be updated again after grads are limited
	petsc_init();

	return;
}

void NavierStokes::solve (int timeStep) {
	initialize_linear_system();
	assemble_linear_system();
	int nIter;
	double rNorm;
	petsc_solve(nIter,rNorm);
	if (Rank==0) cout << "\t" << gid+1 << "\t" << nIter << "\t" << rNorm;
	update_variables();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	calc_limiter();
	mpi_update_ghost_gradients(); // Needs to be updated again after grads are limited
	
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

	qdot.cellStore=false;  qdot.allocate(gid);
	mdot.cellStore=false; mdot.faceStore=true; mdot.allocate(gid);	
	
	// This one is for other solvers to find out contribution from either left or right side
	weightL.cellStore=false; weightL.faceStore=true; weightL.allocate(gid);
	
	return;
}

void NavierStokes::apply_initial_conditions (void) {
	// Loop through each initial condition region and apply sequentially
	int count=input.section("grid",gid).subsection("IC",0).count;
	for (int ic=0;ic<count;++ic) {
		// Store the reference to current IC region
		Subsection &region=input.section("grid",gid).subsection("IC",ic);
		Vec3D regionV=region.get_Vec3D("V");
		double regionRho,regionP,regionT;
		// Assign specified values
		regionP=region.get_double("p");
		if (region.get_double("rho").is_found) { // Rho is specified in the input file
			//  Check if Temperature is also specified
			if (region.get_double("T").is_found) {
				cerr << "Both rho and T can't be specified in initial condition IC_" << ic+1 << " of grid " << gid+1 << endl;
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
	}
	
	// DEBUG
	//for (int c=0;c<grid[gid].cellCount;++c) {
	//	rho.cell(c)=grid[gid].cell[c].centroid[0];
	//	V.cell(c)[0]=grid[gid].cell[c].centroid[0];
	//}
	
	// initialize updates and limiter to zero
	for (int c=0;c<grid[gid].cellCount;++c) {
		for (int i=0;i<5;++i) {
			update[i].cell(c)=0.;
			limiter[i].cell(c)=0.;
		}
	}
	for (int g=0;g<grid[gid].ghostCount;++g) {
		for (int i=0;i<5;++i) {
			update[i].ghost(g)=0.;
			limiter[i].ghost(g)=0.;
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
		gradu.cell(c)=grad[0];
		gradv.cell(c)=grad[1];
		gradw.cell(c)=grad[2];
	 }
	return;
}

void NavierStokes::update_variables(void) {
	
	resP=0.; resV=0.; resT=0.;
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
		
		resP+=update[0].cell(c)*update[0].cell(c);
		resV+=update[1].cell(c)*update[1].cell(c)+update[2].cell(c)*update[2].cell(c)+update[3].cell(c)*update[3].cell(c);
		resT+=update[4].cell(c)*update[4].cell(c);
		
	} // cell loop
	double residuals[3],totalResiduals[3];
	residuals[0]=resP; residuals[1]=resV; residuals[2]=resT;
	
	MPI_Allreduce(&residuals,&totalResiduals,3, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	resP=totalResiduals[0]; resV=totalResiduals[1]; resT=totalResiduals[2];
	
	resP=sqrt(resP); resV=sqrt(resV); resT=sqrt(resT);
	
	// TODO: fix normalization
	/*
	resP/=resP_norm*double(grid[gid].globalCellCount);
	resV/=resV_norm*double(grid[gid].globalCellCount);
	resT/=resT_norm*double(grid[gid].globalCellCount);
	*/
	/*
	if (ps_timeStep==1) {
		resP_first=resP;
		resV_first=resV;
		resT_first=resT;
	}
	 */
	
	return;
}



