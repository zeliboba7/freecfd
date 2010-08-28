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

void NavierStokes::set_bcs(void) {
	
	// Loop through each boundary condition region and apply sequentially
	int count=input.section("grid",gid).subsection("BC",0).count;
	
	for (int b=0;b<count;++b) {
		// Store the reference to current BC region
		Subsection &region=input.section("grid",gid).subsection("BC",b);
		string type=region.get_string("type");
		string kind=region.get_string("kind");
		
		if (kind=="none") bc[gid][b].kind=NONE;
		if (kind=="noreverse") bc[gid][b].kind=NO_REVERSE;
		if (kind=="dampreverse") bc[gid][b].kind=DAMP_REVERSE;
		if (kind=="slip") bc[gid][b].kind=SLIP;
		
		if (region.get_double("p_total").is_found && region.get_double("T_total").is_found) {
			bc[gid][b].kind=STAGNATION;
			p_total.fixedonBC[b]=true;
			p_total.bcValue[b].resize(1);
			p_total.bc(b)=region.get_double("p_total");
			T_total.fixedonBC[b]=true;
			T_total.bcValue[b].resize(1);
			T_total.bc(b)=region.get_double("T_total");
		}
		
		if (region.get_double("p").is_found) {
			p.fixedonBC[b]=true;
			// TODO: If distribution, do the proper thing here
			// p.bc(b,f)=blah;
			// Currently assuming a single uniform value is specified
			// Fix this for all the following once the distributed BC option is implemented
			p.bcValue[b].resize(1);
			p.bc(b)=region.get_double("p");
			bc[gid][b].specified=BC_P;			
			if (region.get_double("T").is_found) {
				if (region.get_double("rho").is_found) {
					cerr << "[E] Thermodynamic state is overspecified for grid " << gid+1 << " boundary condition BC_" << b+1 << endl;
					cerr << "[E] You can only specify up to two of rho, T and p" << endl;
					exit(1);
				}
				// TODO: see above
				T.fixedonBC[b]=true; T.bcValue[b].resize(1);
				T.bc(b)=region.get_double("T");
				// Now get rho from equation of state
				rho.fixedonBC[b]=true; rho.bcValue[b].resize(1);
				rho.bc(b)=material.rho(p.bc(b),T.bc(b));
				bc[gid][b].thermalType=FIXED_T;
				bc[gid][b].specified=BC_STATE;
			} else if (region.get_double("rho").is_found) {
				rho.fixedonBC[b]=true; rho.bcValue[b].resize(1);
				rho.bc(b)=region.get_double("rho");			
				T.fixedonBC[b]=true; T.bcValue[b].resize(1);
				T.bc(b)=material.T(p.bc(b),rho.bc(b));
				bc[gid][b].thermalType=FIXED_T;
				bc[gid][b].specified=BC_STATE;
			} else {
				if (bc[gid][b].type==WALL) bc[gid][b].thermalType=ADIABATIC;
			}
			
		} else if (region.get_double("T").is_found) {
			T.fixedonBC[b]=true; T.bcValue[b].resize(1);
			T.bc(b)=region.get_double("T");
			bc[gid][b].thermalType=FIXED_T;
			bc[gid][b].specified=BC_T;
			if (region.get_double("rho").is_found) {
				rho.fixedonBC[b]=true; rho.bcValue[b].resize(1);
				rho.bc(b)=region.get_double("rho");
				p.fixedonBC[b]=true; p.bcValue[b].resize(1);
				p.bc(b)=material.p(rho.bc(b),T.bc(b));
				bc[gid][b].specified=BC_STATE;
			}
		} else if (region.get_double("rho").is_found) {
			rho.fixedonBC[b]=true; rho.bcValue[b].resize(1);
			rho.bc(b)=region.get_double("rho");
			bc[gid][b].specified=BC_RHO;
		} else if (region.get_double("qdot").is_found) {
			if (bc[gid][b].type==WALL) {
				bc[gid][b].thermalType=FIXED_Q;
				qdot.fixedonBC[b]=true; qdot.bcValue[b].resize(1);
				qdot.bc(b)=region.get_double("qdot");
			}
		} else {
			if (bc[gid][b].type==WALL) {
				bc[gid][b].thermalType=ADIABATIC;
			}
		}
		
		if (bc[gid][b].type==WALL) {
			if (bc[gid][b].kind!=SLIP) {
				V.fixedonBC[b]=true; V.bcValue[b].resize(1); V.bc(b)=0.;
			}
		}
		
		if (bc[gid][b].type==SYMMETRY) bc[gid][b].thermalType=ADIABATIC;
		
		if (region.get_Vec3D("V").is_found) {
			V.fixedonBC[b]=true; V.bcValue[b].resize(1);
			V.bc(b)=region.get_Vec3D("V");
			bc[gid][b].kind=VELOCITY;
		}
		if (region.get_double("mdot").is_found) {
			bc[gid][b].kind=MDOT;
			mdot.fixedonBC[b]=true; mdot.bcValue[b].resize(1);
			mdot.bc(b)=region.get_double("mdot");
		}
				
	}
	
	// TODO: The following is done just to expand the arrays. 
	// If in the future, the bcValue arrays are by default full size, get rid of this
	// Loop through the interfaces
	int b;
	for (int g=0;g<grid.size();++g) {
		for (int i=0;i<interface[g].size();++i) {
			if (interface[g][i].donor_grid==gid) { // If this is a donor
				b=interface[g][i].donor_bc;
				if(interface[g][i].donor_var=="qdot") {
					qdot.bcValue[b].resize(grid[gid].boundaryFaceCount[b][grid[gid].Rank]);
				}
				// Flowfield variables do not need to be allocated at the bc's
			} // if donor
			else if (interface[g][i].recv_grid==gid) { // If this is a receiver
				b=interface[g][i].recv_bc;
				if(interface[g][i].recv_var=="qdot") {
					qdot.bcValue[b].resize(grid[gid].boundaryFaceCount[b][grid[gid].Rank]);
					for (int j=1;j<qdot.bcValue[b].size();++j) qdot.bcValue[b][j]=qdot.bcValue[b][0];
				}
				if(interface[g][i].recv_var=="T") {
					T.bcValue[b].resize(grid[gid].boundaryFaceCount[b][grid[gid].Rank]);
					for (int j=1;j<T.bcValue[b].size();++j) T.bcValue[b][j]=T.bcValue[b][0];
				}
			} // if receiver
		}
	}

	return;
}
