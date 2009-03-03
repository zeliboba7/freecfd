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
#include "flamelet.h"

void Flamelet_Table::read(string fileName) {
	//Open file
	// Read nZ,nZvar,nChi
	nZ=1;nZvar=1;nChi=1;
	// Allocate table containers
	rho.resize(nZ);
	for (int i=0;i<nZ;++i) {
		rho[i].resize(nZvar);
		for (int j=0;j<nZvar;++j) {
			rho[i][j].resize(nChi);
		}
	}	
	// Copy for other variables
	T=rho;
	viscosity=rho;
	diffusivity=rho;
	conductivity=rho;
	
	return; 
} // end Flamelet_Table::read

double Flamelet_Table::get_rho(double &Z, double &Zvar, double &Chi) {
	return 2.5e2*pow(Z,4)-5.4e2*pow(Z,3)+3.9e2*pow(Z,2)-1.1e2*Z+13.0;
}

double Flamelet_Table::get_temperature(double &Z, double &Zvar, double &Chi) {
	return -4.65e4*pow(Z,4)+1.12e5*pow(Z,3)-9.26e4*pow(Z,2)+2.67e4*Z+1.03e3;
}