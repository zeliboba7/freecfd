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
#include "fstream"

void Flamelet_Table::read(string fileName) {
	string tag;
	int nVar;

	if (Rank==0) cout << "[I] Reading the flamelet table" << endl;
	
	//Open file
	ifstream chemTable;
	chemTable.open(fileName.c_str(), ios::in);
	if (!chemTable.is_open()) {
		if (Rank==0) cerr << "[E] Flamelet table file " << fileName  << " couldn't be opened" << endl;
		exit(1);
	}

	chemTable >> nZ; if (Rank==0) cout << "[I] ...Number of Z points: " << nZ << endl;
	chemTable >> nZvar; if (Rank==0) cout << "[I] ...Number of Zvar points: " << nZvar << endl;
	chemTable >> nChi; if (Rank==0) cout << "[I] ...Number of Chi points: " << nChi << endl;
	chemTable >> nVar;
	
	Z.resize(nZ);
	Zvar.resize(nZvar);
	Chi.resize(nChi);
	weights.resize(6);

	rho.resize(nZ);
	T.resize(nZ);
	viscosity.resize(nZ);
	for (int i=0;i<nZ;++i) {
		rho[i].resize(nZvar);
		T[i].resize(nZvar);
		viscosity[i].resize(nZvar);
		for (int j=0;j<nZvar;++j) {
			rho[i][j].resize(nChi);
			T[i][j].resize(nChi);
			viscosity[i][j].resize(nChi);
		}
	}	
	for (int i=0;i<nZ;++i) chemTable >> Z[i];
	for (int i=0;i<nZvar;++i) chemTable >> Zvar[i];
	for (int i=0;i<nChi;++i) chemTable >> Chi[i];
	
	double dummy;
	
	for (int i=0;i<nVar;++i){
		chemTable >> tag;
		if (Rank==0) cout << "[I] ...Found variable: " << tag << endl;
		if (tag=="RHO") for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> rho[j][k][m]; 
		else if (tag=="T")for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> T[j][k][m]; 
		else if (tag=="VISC")for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> viscosity[j][k][m]; 
		else for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> dummy;	
      	}
	
	return; 
} // end Flamelet_Table::read

void Flamelet_Table::get_weights(double &Z_in, double &Zvar_in, double &Chi_in) {
	int i1,i2,i3;
	i1=-1;
	if (Z_in< Z[0]) {
		i1=0;
		weights[0]=1.0;
	} else if(Z_in>=Z[1]) {
		i1=nZ-1;
		weights[0]=0.0;
	}else {for (int i=0;i<nZ-1;++i) if (Z_in <=Z[i+1]) { i1=i; break;}
	weights[0]=1.0-(Z_in-Z[i1])/(Z[i1+1]-Z[i1]);	
	}
	weights[1]=1.0-weights[0];
	//allow for non-equdistant grid spacing for x2=Zvar
	i2=-1;
	if (Zvar_in< Zvar[0]) {
		i2=0;
		weights[2]=1.0;
	} else if(Zvar_in>=Zvar[1]) {
		i2=nZvar-1;
		weights[2]=0.0;
	}else {for (int i=0;i<nZvar-1;++i) if (Zvar_in <=Zvar[i2+1]) { i2=i; break;}
		weights[2]=1.0-(Zvar_in-Zvar[i2])/(Zvar[i2+1]-Zvar[i2]);	
	}
	weights[3]=1.0-weights[2];
	//allow for non-equdistant grid spacing for x3=chiMean
	i3=-1;
	if (Chi_in< Chi[0]) {
		i3=0;
		weights[4]=1.0;
	} else if(Chi_in>Chi[0]) {
		i3=nChi-1;
		weights[4]=0.0;
	}else {for (int i=0;i<nChi-1;++i) if (Chi_in <=Chi[i+1]) { i3=i; break;}
		weights[4]=1.0-(Chi_in-Chi[i3])/(Chi[i3+1]-Chi[i3]);	
	}
	weights[5]=1.0-weights[4];
	// A little check
	if ((i1+1)*(i2+1)*(i3+1)==0){
		cout<<"Error in Lookup table"<<endl;
		exit(1);
	}
	
	
}


double Flamelet_Table::get_rho(double &Z, double &Zvar, double &Chi,bool refreshWeights) {
	if (refreshWeights) get_weights(Z,Zvar,Chi);
	return  weights[4]*( weights[2]*( weights[0]*rho[i1][i2][i3]
			+weights[1]*rho[i1+1][i2][i3] )
			+weights[3]*( weights[0]*rho[i1][i2+1][i3]
			+weights[1]*rho[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*rho[i1][i2][i3+1]   
			+weights[1]*rho[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*rho[i1][i2+1][i3+1] 
			+weights[1]*rho[i1+1][i2+1][i3+1] ) );
	
	//return 2.5e2*pow(Z,4)-5.4e2*pow(Z,3)+3.9e2*pow(Z,2)-1.1e2*Z+13.0;
}

double Flamelet_Table::get_temperature(double &Z, double &Zvar, double &Chi,bool refreshWeights) {
	if (refreshWeights) get_weights(Z,Zvar,Chi);
	return  weights[4]*( weights[2]*( weights[0]*T[i1][i2][i3]
			+weights[1]*T[i1+1][i2][i3] )
			+weights[3]*( weights[0]*T[i1][i2+1][i3]
			+weights[1]*T[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*T[i1][i2][i3+1]   
			+weights[1]*T[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*T[i1][i2+1][i3+1] 
			+weights[1]*T[i1+1][i2+1][i3+1] ) );
	
	//return -4.65e4*pow(Z,4)+1.12e5*pow(Z,3)-9.26e4*pow(Z,2)+2.67e4*Z+1.03e3;
}