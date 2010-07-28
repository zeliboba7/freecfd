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
#include "rans.h"

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

void RANS::initialize (void) {
	
	nVars=2;
	rtol=input.section("grid",gid).subsection("turbulence").get_double("relativetolerance");
	abstol=input.section("grid",gid).subsection("turbulence").get_double("absolutetolerance");
	maxits=input.section("grid",gid).subsection("turbulence").get_int("maximumiterations");
	
	if (input.section("grid",gid).subsection("turbulence").get_string("model")=="k-omega") {
		model=KOMEGA;
	} else if (input.section("grid",gid).subsection("turbulence").get_string("model")=="k-epsilon") {
		model=KEPSILON;
	} else if (input.section("grid",gid).subsection("turbulence").get_string("model")=="sst") {
		model=SST;
	}

	if (input.section("grid",gid).subsection("turbulence").get_string("order")=="first") {
		order=FIRST;
	} else if (input.section("grid",gid).subsection("turbulence").get_string("order")=="second") {
		order=SECOND;
	}
	
	kLowLimit=input.section("grid",0).subsection("turbulence").get_double("klowlimit");
	kHighLimit=input.section("grid",0).subsection("turbulence").get_double("khighlimit");
	omegaLowLimit=input.section("grid",0).subsection("turbulence").get_double("omegalowlimit");
	viscosityRatioLimit=input.section("grid",0).subsection("turbulence").get_double("viscosityratiolimit");
	Pr_t=input.section("grid",0).subsection("turbulence").get_double("turbulentPr");

	
	mpi_init();
	material.set(gid);
	create_vars();
	apply_initial_conditions();
	set_bcs();
	update_eddy_viscosity();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	petsc_init();

	
	return;
}

void RANS::solve (int timeStep) {

	terms();
	int nIter;
	double rNorm;
	petsc_solve(nIter,rNorm);
	if (Rank==0) cout << "\t" << gid+1 << "\t" << nIter << "\t" << rNorm << "\t";
	update_variables();
	update_eddy_viscosity();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();

	return;
}

void RANS::create_vars (void) {
	// Allocate variables
	// Default option is to store on cell centers and ghosts only
	k.allocate(gid);
	omega.allocate(gid);
	mu_t.allocate(gid);
	strainRate.allocate(gid);
	gradk.allocate(gid);
	gradomega.allocate(gid);
	update.resize(2);
	for (int i=0; i<2; ++i) update[i].allocate(gid);
	
	return;
}

void RANS::apply_initial_conditions (void) {
	// Loop through each initial condition region and apply sequentially
	int count=input.section("grid",gid).subsection("IC",0).count;
	for (int ic=0;ic<count;++ic) {
		// Store the reference to current IC region
		Subsection &region=input.section("grid",gid).subsection("IC",ic);
		Vec3D regionV=region.get_Vec3D("V");
		double intensity,viscRatio;
		// Assign specified values
		intensity=region.get_double("turbulenceintensity");
		viscRatio=region.get_double("eddyviscosityratio");

		// If region is specified with a box method
		if (region.get_string("region")=="box") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the box region
				if (withinBox(grid[gid].cell[c].centroid,region.get_Vec3D("corner_1"),region.get_Vec3D("corner_2"))) {
					k.cell(c)=intensity*fabs(ns[gid].V.cell(c));
					k.cell(c)*=1.5*k.cell(c);
					omega.cell(c)=viscRatio*material.viscosity(ns[gid].T.cell(c))/(ns[gid].rho.cell(c)*k.cell(c));
				}
			}
		} else if (region.get_string("region")=="cylinder") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the cylinder region
				Vec3D axisDirection=region.get_Vec3D("axisdirection");
				axisDirection=axisDirection.norm();
				if (withinCylinder(grid[gid].cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"),axisDirection,region.get_double("height"))) {
					k.cell(c)=intensity*fabs(ns[gid].V.cell(c));
					k.cell(c)*=1.5*k.cell(c);
					omega.cell(c)=viscRatio*material.viscosity(ns[gid].T.cell(c))/(ns[gid].rho.cell(c)*k.cell(c));
				}
			}
		} else if (region.get_string("region")=="sphere") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the sphere region
				if (withinSphere(grid[gid].cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"))) {
					k.cell(c)=intensity*fabs(ns[gid].V.cell(c));
					k.cell(c)*=1.5*k.cell(c);
					omega.cell(c)=viscRatio*material.viscosity(ns[gid].T.cell(c))/(ns[gid].rho.cell(c)*k.cell(c));
				}
			}
		}
	}

	// initialize updates and limiter
	for (int c=0;c<grid[gid].cellCount;++c) for (int i=0;i<2;++i) update[i].cell(c)=0.;

	for (int g=0;g<grid[gid].ghostCount;++g) for (int i=0;i<2;++i) update[i].ghost(g)=0.;

	return;
}

void RANS::calc_cell_grads (void) {
	for (int c=0;c<grid[gid].cellCount;++c) {
		gradk.cell(c)=k.cell_gradient(c);
		gradomega.cell(c)=omega.cell_gradient(c);
	}
	return;
}

void RANS::update_variables(void) {
	
	//resK=0.; resOmega=0.;
	int counter=0;
	double mu,new_mu_t;
	for (int c=0;c<grid[gid].cellCount;++c) {
		for (int i=0;i<2;++i) {
			if (isnan(update[i].cell(c)) || isinf(update[i].cell(c))) {
				cerr << "[E] Divergence detected!...exiting" << endl;
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
		// Limit the update so that k and omega doesn't end up out of limits
		update[0].cell(c)=max(-1.*(k.cell(c)-kLowLimit),update[0].cell(c));
		update[0].cell(c)=min((kHighLimit-k.cell(c)),update[0].cell(c));
		
		update[1].cell(c)=max(-1.*(omega.cell(c)-omegaLowLimit),update[1].cell(c));
		
		new_mu_t=ns[gid].rho.cell(c)*(k.cell(c)+update[0].cell(c))/(omega.cell(c)+update[1].cell(c));
		mu=material.viscosity(ns[gid].T.cell(c));
		if (new_mu_t/mu>viscosityRatioLimit) {
			counter++; 
			double under_relax;
			double limit_nu=viscosityRatioLimit*mu/ns[gid].rho.cell(c);
			under_relax=(limit_nu*omega.cell(c)-k.cell(c))/(update[0].cell(c)-limit_nu*update[1].cell(c)+1.E-8);
			under_relax=0.9*max(1.,under_relax);
			update[0].cell(c)*=under_relax;
			update[1].cell(c)*=under_relax;
		}
		
		k.cell(c) += update[0].cell(c);
		omega.cell(c)+= update[1].cell(c);
		
		//resK+=update[0].cell(c)*update[0].cell(c);
		//resOmega+=update[1].cell(c)*update[1].cell(c);
		
	} // cell loop
	if (counter>0) cout << "[I] Update of k and omega is limited due to viscosityRatioLimit constraint for " << counter << " cells" << endl;

	/*
	double residuals[2],totalResiduals[2];
	residuals[0]=resK; residuals[1]=resOmega;
	if (np!=1) {
		MPI_Reduce(&residuals,&totalResiduals,2, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		resK=totalResiduals[0]; resOmega=totalResiduals[1];
	}
	resK=sqrt(resK); resOmega=sqrt(resOmega);
	resK/=resK_norm*double(grid.globalCellCount);
	resOmega/=resOmega_norm*double(grid.globalCellCount);
	*/
	return;
}
