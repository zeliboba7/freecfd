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
extern RANS_Model komega,kepsilon;

void RANS::update_eddy_viscosity(void) {
	
	double arg1,arg2,arg3,F2;
	double mu;
	double a1=0.31; // SST a1 value
	double turbulent_length_scale;

	for (int c=0;c<grid[gid].cell.size();++c) {
		mu=ns[gid].material.viscosity(ns[gid].T.cell(c));
		if (model==SST) {
			arg1=2.*sqrt(k.cell(c)+1.e-15)/(kepsilon.beta_star*omega.cell(c)*grid[gid].cell[c].closest_wall_distance);
			arg2=500.*mu/(ns[gid].rho.cell(c)*omega.cell(c)
					*grid[gid].cell[c].closest_wall_distance*grid[gid].cell[c].closest_wall_distance);
			arg3=max(arg1,arg2);
			F2=tanh(arg3*arg3);
			mu_t.cell(c)=a1*ns[gid].rho.cell(c)*k.cell(c)/max(a1*omega.cell(c),strainRate.cell(c)*F2);
			mu_t.cell(c)=min(viscosityRatioLimit*mu,mu_t.cell(c));
			//mu_t.cell(c)=ns[gid].rho.cell(c)*k.cell(c)/omega.cell(c);
			//turbulent_length_scale=sqrt((k.cell(c)+1.e-15))/max(omega.cell(c),strainRate.cell(c)*F2/a1);
			//turbulent_length_scale/=0.083;
		} else {
			mu_t.cell(c)=ns[gid].rho.cell(c)*k.cell(c)/omega.cell(c);
			mu_t.cell(c)=min(viscosityRatioLimit*mu,mu_t.cell(c));
			//turbulent_length_scale=sqrt(k.cell(c)+1.e-15)/omega.cell(c);
			//turbulent_length_scale/=0.083;
		}
	}
	
	return;
	
} 



