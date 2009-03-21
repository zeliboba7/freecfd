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
#include <iomanip>


Species::Species(string name_in) {
	if (name_in=="O2") Mw=31.9988*1.e-3;
	else if (name_in=="H2") Mw=2.01594*1.e-3;
	else if (name_in=="H2O") Mw=18.01534*1.e-3;
	else if (name_in=="OH") Mw=17.00737*1.e-3;
	else {
	 	if (Rank==0) cerr << "[E] Species " << name_in << " is not recognized" << endl;
		exit(1);
	}
	name=name_in;
	return;
}

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
	mu.resize(nZ);
	diffusivity.resize(nZ);
	Y.resize(nZ);
	for (int i=0;i<nZ;++i) {
		rho[i].resize(nZvar);
		T[i].resize(nZvar);
		mu[i].resize(nZvar);
		diffusivity[i].resize(nZvar);
		Y[i].resize(nZvar);
		for (int j=0;j<nZvar;++j) {
			rho[i][j].resize(nChi);
			T[i][j].resize(nChi);
			mu[i][j].resize(nChi);
			diffusivity[i][j].resize(nChi);
			Y[i][j].resize(nChi);
		}
	}	
	for (int i=0;i<nZ;++i) chemTable >> Z[i];
	for (int i=0;i<nZvar;++i) chemTable >> Zvar[i];
	for (int i=0;i<nChi;++i) chemTable >> Chi[i];
	
	double dummy;
	
	for (int i=0;i<nVar;++i){
		chemTable >> tag;
		if (Rank==0) cout << "[I] ...Found variable: " << tag << endl;
		if (tag=="RHO") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> rho[j][k][m];
		} else if (tag=="T") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> T[j][k][m]; 
		} else if (tag=="VISC") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> mu[j][k][m]; 
		} else if (tag=="DIFF") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> diffusivity[j][k][m]; 
		} else if (tag.substr(0,2)=="Y_") {
			Species temp(tag.substr(2));
			species.push_back(temp);
			for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) {
				Y[j][k][m].resize(species.size());
				chemTable >> Y[j][k][m][species.size()-1];
			}
		} else {
			for (int j=0;j<nZ;++j) for (int k=0;k<nZvar;++k) for (int m=0;m<nChi;++m) chemTable >> dummy;
		}		
      	}
	
	return; 
} // end Flamelet_Table::read

void Flamelet_Table::get_weights(double &Z_in, double &Zvar_in, double &Chi_in) {
	
	i1=-1;
	if (Z_in<=Z[0]) {
		i1=0;
		weights[0]=1.0;
	} else if (Z_in>=Z[nZ-1]) {
		i1=nZ-2;
		weights[0]=0.0;
	} else {
		for (int i=0;i<nZ-1;++i) {
			if (Z_in<=Z[i+1]) { 
				i1=i; break;
			}
		}
		weights[0]=1.0-(Z_in-Z[i1])/(Z[i1+1]-Z[i1]);	
	}
	weights[1]=1.0-weights[0];
	//allow for non-equdistant grid spacing for x2=Zvar
	i2=-1;
	if (Zvar_in<=Zvar[0]) {
		i2=0;
		weights[2]=1.0;
	} else if (Zvar_in>=Zvar[nZvar-1]) {
		i2=nZvar-2;
		weights[2]=0.0;
	} else {
		for (int i=0;i<nZvar-1;++i) {
			if (Zvar_in<=Zvar[i+1]) { 
				i2=i; break;
			}
		}
		weights[2]=1.0-(Zvar_in-Zvar[i2])/(Zvar[i2+1]-Zvar[i2]);	
	}
	weights[3]=1.0-weights[2];
	//allow for non-equdistant grid spacing for x3=chiMean
	i3=-1;
	if (Chi_in<=Chi[0]) {
		i3=0;
		weights[4]=1.0;
	} else if (Chi_in>=Chi[nChi-1]) {
		i3=nChi-2;
		weights[4]=0.0;
	} else {
		for (int i=0;i<nChi-1;++i) {
			if (Chi_in<=Chi[i+1]) {
				i3=i; break;
			}
		}
		weights[4]=1.0-(Chi_in-Chi[i3])/(Chi[i3+1]-Chi[i3]);	
	}
	weights[5]=1.0-weights[4];
	// A little check
	if ((i1+1)*(i2+1)*(i3+1)==0){
		cout<<"Error in Lookup table"<<endl;
		exit(1);
	}
	
	return;
}


double Flamelet_Table::get_rho(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return weights[4]*( weights[2]*( weights[0]*rho[i1][i2][i3]
			+weights[1]*rho[i1+1][i2][i3] )
			+weights[3]*( weights[0]*rho[i1][i2+1][i3]
			+weights[1]*rho[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*rho[i1][i2][i3+1]   
			+weights[1]*rho[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*rho[i1][i2+1][i3+1] 
			+weights[1]*rho[i1+1][i2+1][i3+1] ) );

}

double Flamelet_Table::get_T(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return weights[4]*( weights[2]*( weights[0]*T[i1][i2][i3]
			+weights[1]*T[i1+1][i2][i3] )
			+weights[3]*( weights[0]*T[i1][i2+1][i3]
			+weights[1]*T[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*T[i1][i2][i3+1]   
			+weights[1]*T[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*T[i1][i2+1][i3+1] 
			+weights[1]*T[i1+1][i2+1][i3+1] ) )-Tref;
}

void Flamelet_Table::get_rho_T_comp(double &p_in, double &Z_in, double &Zvar_in, double &Chi_in,double &rho_out, double &T_out) {
	
	double rho_table,T_table,p_table,Mw_table;
	
	T_table=get_T(Z_in,Zvar_in,Chi_in);
	rho_table=get_rho(Z_in,Zvar_in,Chi_in,false);
	Mw_table=get_Mw(Z_in,Zvar_in,Chi_in,false);
	p_table=rho_table*UNIV_GAS_CONST/Mw_table*T_table;
	
	// Assume isentropic compression or expansion from p_table to p_in and correct temperature and density
	rho_out=rho_table*pow((p_in+Pref)/p_table,1./Gamma);
	T_out=(p_in+Pref)*Mw_table/(UNIV_GAS_CONST*rho_out)-Tref;
	
	//Deactivate
// 	rho_out=rho_table;
// 	T_out=T_table-Tref;
	//T_out=pow( pow(p_table,Gamma-1.)*pow(T_table,-Gamma)/pow(p_in+Pref,Gamma-1.) , -1./Gamma);
	//rho_out=(p_in+Pref)*Mw_table/(UNIV_GAS_CONST*T_out);
	//T_out-=Tref;
	return;
}

double Flamelet_Table::get_mu(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return  weights[4]*( weights[2]*( weights[0]*mu[i1][i2][i3]
			+weights[1]*mu[i1+1][i2][i3] )
			+weights[3]*( weights[0]*mu[i1][i2+1][i3]
			+weights[1]*mu[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*mu[i1][i2][i3+1]   
			+weights[1]*mu[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*mu[i1][i2+1][i3+1] 
			+weights[1]*mu[i1+1][i2+1][i3+1] ) );
}

double Flamelet_Table::get_diffusivity(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return  weights[4]*( weights[2]*( weights[0]*diffusivity[i1][i2][i3]
			+weights[1]*diffusivity[i1+1][i2][i3] )
			+weights[3]*( weights[0]*diffusivity[i1][i2+1][i3]
			+weights[1]*diffusivity[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*diffusivity[i1][i2][i3+1]   
			+weights[1]*diffusivity[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*diffusivity[i1][i2+1][i3+1] 
			+weights[1]*diffusivity[i1+1][i2+1][i3+1] ) );
}

double Flamelet_Table::get_Mw(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	double Mw=0.;
	double Mw_s,Y_s;
	double Y_sum=0.;
	for (int s=0;s<species.size();++s) {
		Y_s=weights[4]*( weights[2]*( weights[0]*Y[i1][i2][i3][s]
				+weights[1]*Y[i1+1][i2][i3][s] )
				+weights[3]*( weights[0]*Y[i1][i2+1][i3][s]
				+weights[1]*Y[i1+1][i2+1][i3][s] ) ) 
				+weights[5]*( weights[2]*( weights[0]*Y[i1][i2][i3+1][s] 
				+weights[1]*Y[i1+1][i2][i3+1][s] ) 
				+weights[3]*( weights[0]*Y[i1][i2+1][i3+1][s]
				+weights[1]*Y[i1+1][i2+1][i3+1][s] ) );
		Mw+=Y_s/species[s].Mw;
		Y_sum+=Y_s;
	}
	Mw=1./Mw;
	Mw/=Y_sum; // small fix just in case sum of all the species mass fractions are not exactly 1
	
	return Mw;
}

double Flamelet_Table::get_drho_dZ(double &Z_in, double &Zvar_in, double &Chi_in) {
	
	get_weights(Z_in,Zvar_in,Chi_in);
	// Know we know i1
	double Z_minus,Z_plus;
	if (i1<Z.size()) {
		Z_minus=Z[i1];
		Z_plus=Z[i1+1];
	} else {
		Z_minus=Z[i1-1];
	 	Z_plus=Z[i1];
	}
	
	return (get_rho(Z_plus,Zvar_in,Chi_in)-get_rho(Z_minus,Zvar_in,Chi_in))/(Z_plus-Z_minus);
}



double Flamelet_Table::get_drho_dZvar(double &Z_in, double &Zvar_in, double &Chi_in) {
	
	get_weights(Z_in,Zvar_in,Chi_in);
	// Know we know i2
	double Zvar_minus,Zvar_plus;
	if (i2<Zvar.size()) {
		Zvar_minus=Z[i2];
		Zvar_plus=Z[i2+1];
	} else {
		Zvar_minus=Z[i2-1];
		Zvar_plus=Z[i2];
	}
	return (get_rho(Z_in,Zvar_plus,Chi_in)-get_rho(Z_in,Zvar_minus,Chi_in))/(Zvar_plus-Zvar_minus);
}