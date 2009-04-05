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
	chemTable >> nS; if (Rank==0) cout << "[I] ...Number of S points: " << nS << endl;
	chemTable >> nChi; if (Rank==0) cout << "[I] ...Number of Chi points: " << nChi << endl;
	chemTable >> nVar; if (Rank==0) cout << "[I] ...Number of variables: " << nVar << endl;
	
	Z.resize(nZ);
	S.resize(nS);
	Chi.resize(nChi);
	weights.resize(6);

	rho.resize(nZ);
	T.resize(nZ);
	mu.resize(nZ);
	diffusivity.resize(nZ);
	c_p.resize(nZ);
	Mw.resize(nZ);
	//Y.resize(nZ);
	for (int i=0;i<nZ;++i) {
		rho[i].resize(nS);
		T[i].resize(nS);
		mu[i].resize(nS);
		diffusivity[i].resize(nS);
		c_p[i].resize(nS);
		Mw[i].resize(nS);
		//Y[i].resize(nS);
		for (int j=0;j<nS;++j) {
			rho[i][j].resize(nChi);
			T[i][j].resize(nChi);
			mu[i][j].resize(nChi);
			diffusivity[i][j].resize(nChi);
			c_p[i][j].resize(nChi);
			Mw[i][j].resize(nChi);
			//Y[i][j].resize(nChi);
		}
	}	
	for (int i=0;i<nZ;++i) chemTable >> Z[i];
	for (int i=0;i<nS;++i) chemTable >> S[i];
	for (int i=0;i<nChi;++i) chemTable >> Chi[i];
	
	double dummy;
	
	for (int i=0;i<nVar;++i){
		chemTable >> tag;
		if (Rank==0) cout << "[I] ...Found variable: " << tag << endl;
		if (tag=="RHO") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> rho[j][k][m];
		} else if (tag=="T") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> T[j][k][m]; 
		} else if (tag=="VISC") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> mu[j][k][m]; 
		} else if (tag=="DIFF") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> diffusivity[j][k][m]; 
// 		} else if (tag.substr(0,2)=="Y_") {
// 			Species temp(tag.substr(2));
// 			species.push_back(temp);
// 			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) {
// 				Y[j][k][m].resize(species.size());
// 				chemTable >> Y[j][k][m][species.size()-1];
// 			}
		} else if (tag=="C_P") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> c_p[j][k][m];
		} else if (tag=="MW") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) {
				chemTable >> Mw[j][k][m];
				Mw[j][k][m]/=1000.;
			}
		} else {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> dummy;
		}		
      	}
	
	return; 
} // end Flamelet_Table::read

void Flamelet_Table::get_weights(double &Z_in, double &Zvar_in, double &Chi_in) {
	
	// Convert Zvar_in to S_in
	double S_in;
	
	if ( fabs(Z_in)<1.e-8 || fabs(Z_in-1.)<1.e-8) {
		S_in=0.;
	} else {
		S_in=Zvar_in/(Z_in*(1.-Z_in));
	}
	
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
	
	i2=-1;
	if (S_in<=S[0]) {
		i2=0;
		weights[2]=1.0;
	} else if (S_in>=S[nS-1]) {
		i2=nS-2;
		weights[2]=0.0;
	} else {
		for (int i=0;i<nS-1;++i) {
			if (S_in<=S[i+1]) { 
				i2=i; break;
			}
		}
		weights[2]=1.0-(S_in-S[i2])/(S[i2+1]-S[i2]);	
	}
	weights[3]=1.0-weights[2];
	
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
	
	return weights[4]*( weights[2]*( weights[0]*rho[i1][i2][i3]+weights[1]*rho[i1+1][i2][i3] )
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

double Flamelet_Table::get_c_p(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return  weights[4]*( weights[2]*( weights[0]*c_p[i1][i2][i3]
			+weights[1]*c_p[i1+1][i2][i3] )
			+weights[3]*( weights[0]*c_p[i1][i2+1][i3]
			+weights[1]*c_p[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*c_p[i1][i2][i3+1]
			+weights[1]*c_p[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*c_p[i1][i2+1][i3+1] 
			+weights[1]*c_p[i1+1][i2+1][i3+1] ) );
}

double Flamelet_Table::get_Mw(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
// 	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
// 	
// 	return  weights[4]*( weights[2]*( weights[0]*Mw[i1][i2][i3]
// 			+weights[1]*Mw[i1+1][i2][i3] )
// 			+weights[3]*( weights[0]*Mw[i1][i2+1][i3]
// 			+weights[1]*Mw[i1+1][i2+1][i3] ) ) 
// 			+weights[5]*( weights[2]*( weights[0]*Mw[i1][i2][i3+1]   
// 			+weights[1]*Mw[i1+1][i2][i3+1] ) 
// 			+weights[3]*( weights[0]*Mw[i1][i2+1][i3+1] 
// 			+weights[1]*Mw[i1+1][i2+1][i3+1] ) );
// 	
	 return UNIV_GAS_CONST/(Pref/(get_rho(Z_in,Zvar_in,Chi_in,refreshWeights)*get_T(Z_in,Zvar_in,Chi_in,refreshWeights)));
	
	
}

double Flamelet_Table::get_gamma(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	double c_p_table=get_c_p(Z_in,Zvar_in,Chi_in,refreshWeights);
	double R=UNIV_GAS_CONST/get_Mw(Z_in,Zvar_in,Chi_in,false);
		
	return c_p_table/(c_p_table-R);
}

void Flamelet_Table::get_rho_T_comp(double &p_in, double &Z_in, double &Zvar_in, double &Chi_in,double &rho_out, double &T_out) {
	
	double rho_table,T_table,p_table,Mw_table;
	
	T_table=get_T(Z_in,Zvar_in,Chi_in);
	rho_table=get_rho(Z_in,Zvar_in,Chi_in,false);
	Gamma=get_gamma(Z_in,Zvar_in,Chi_in,false);
	Mw_table=get_Mw(Z_in,Zvar_in,Chi_in,false);
	p_table=rho_table*UNIV_GAS_CONST/Mw_table*T_table;
	
	// Assume isentropic compression or expansion from p_table to p_in and correct temperature and density
	rho_out=rho_table*pow((p_in+Pref)/p_table,1./Gamma);
	T_out=(p_in+Pref)*Mw_table/(UNIV_GAS_CONST*rho_out)-Tref;

	//Deactivate
	rho_out=rho_table;
	T_out=T_table-Tref;

	return;
}

double Flamelet_Table::get_drho_dZ(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	// Now we know i1
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

double Flamelet_Table::get_drho_dZvar(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	// Know we know i2
	double Zvar_minus,Zvar_plus;
	double Zvar2S=Z_in*(1.-Z_in);
	if (i2<S.size()) {
		Zvar_minus=S[i2]*Zvar2S;
		Zvar_plus=S[i2+1]*Zvar2S;
	} else {
		Zvar_minus=S[i2-1]*Zvar2S;
		Zvar_plus=S[i2]*Zvar2S;
	}
	if (fabs(Zvar_plus-Zvar_minus)<1.e-8) return 0.;
	return (get_rho(Z_in,Zvar_plus,Chi_in)-get_rho(Z_in,Zvar_minus,Chi_in))/(Zvar_plus-Zvar_minus);
	
}