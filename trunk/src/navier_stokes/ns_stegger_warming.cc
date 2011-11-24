#include "ns.h"

void mtx_out(double L[5][5],int dim);

inline
void mtx_mult(double A1[5][5],double A2[5][5],double R[5][5]){
	for(int i=0;i<5;i++){
	   for(int j=0;j<5;j++){
	      R[i][j]=0.;
	      for(int k=0;k<5;k++)
	         R[i][j]+=A1[i][k]*A2[k][j];
	   }
	}
	return;
}

inline 
void mtx_vec_mult(double A[5][5],double B[5],double R[5][5]) {
	for (int i=0;i<5;++i) for (int j=0;j<5;++j) R[i][j]+=A[i][j]*B[j];
	return;
}

void Stegger_Warming_flux(NS_Cell_State &left,NS_Cell_State &right,double diss_factor,double closest_wall_distance,double wdiss,double bl_height,MATERIAL &material,double fluxNormal[],double &weightL) { 

	double alpha=6.; // Some problems may need larger alpha
	double rho,p,T,a,H,V2,a2;
	Vec3D V;
	double RW=UNIV_GAS_CONST/material.Mw;
	double gamma_s,beta,eta,Cp;
	double w,deltap,weightR,eps;
	double signal;
	double Lambda[5],L[5][5],R[5][5],Q[5],Jacob[5][5];

	for (int i=0;i<5;++i) fluxNormal[i]=0.;

//	if (closest_wall_distance<bl_height) w=0.5;
//	else {
		deltap=fabs(left.p-right.p)/min(left.p,right.p);
		w=0.5/(alpha*alpha*deltap*deltap+1.);
//	}

	for (int h=0;h<=1;h++) {
	signal=pow(-1.,h);

	if (signal>0.) weightL=1.-w;
	else weightL=w;
	weightR=1.-weightL;

	p=weightL*left.p+weightR*right.p;
	V=weightL*left.Vn+weightR*right.Vn;
	V2=fabs(V); V2*=V2;
	T=weightL*left.T+weightR*right.T;
	rho=weightL*left.rho+weightR*right.rho;
//	rho=material.rho(p,T);
	Cp=material.Cp(T); 
	beta=material.gamma*RW/Cp;
	a=sqrt((1.+beta)*p/rho);
//	a=material.a(p,T); 
	a2=a*a;
	H=0.5*V2+a2/(material.gamma-1.);
	gamma_s=RW*T+beta*(0.5*V2-Cp*T/material.gamma);
	eta=(beta*V2-gamma_s)/(beta*a2);

	for (int i=0;i<5;++i) { 
		Lambda[i]=0.;
		for (int j=0;j<5;++j) {
			L[i][j]=0.;
			R[i][j]=0.;
			Jacob[i][i]=0.;
		}
	}

	L[0][0]=1.0/a2;
	L[0][3]=L[0][4]=0.5/a2;

	L[1][0]=V[0]/a2;
	L[2][0]=V[1]/a2;
	L[3][0]=V[2]/a2;

	L[1][3]=0.5*(V[0]+a)/a2;
	L[1][4]=0.5*(V[0]-a)/a2;
	L[2][1]=1.;
	L[2][3]=L[2][4]=0.5*V[1]/a2;
	L[3][2]=1.;
	L[3][3]=L[3][4]=0.5*V[2]/a2;
	L[4][0]=eta;
	L[4][1]=V[1];
	L[4][2]=V[2];
	L[4][3]=0.5*(H+a*V[0])/a2;
	L[4][4]=0.5*(H-a*V[0])/a2;

	if (closest_wall_distance>bl_height)	eps=wdiss*(a+fabs(V));
	else				eps=wdiss*diss_factor*(a+fabs(V[0]));

	Lambda[0]=Lambda[1]=Lambda[2]=0.5*(V[0]+signal*sqrt(V[0]*V[0]+eps*eps));
	Lambda[3]=0.5*(V[0]+a+signal*sqrt((V[0]+a)*(V[0]+a)+eps*eps));
	Lambda[4]=0.5*(V[0]-a+signal*sqrt((V[0]-a)*(V[0]-a)+eps*eps));

	mtx_vec_mult(L,Lambda,R);
	
	for (int i=0;i<5;++i) for (int j=0;j<5;++j) L[i][j]=R[i][j];

	for (int i=0;i<5;++i) for (int j=0;j<5;++j) R[i][j]=0.;

	R[0][0]=a2-gamma_s;
	R[0][1]=beta*V[0];
	R[0][2]=beta*V[1];
	R[0][3]=beta*V[2];
	R[0][4]=-beta;
	R[1][0]=-V[1];
	R[1][2]=1.;
	R[2][0]=-V[2];
	R[2][3]=1.;
	R[3][0]=gamma_s-V[0]*a;
	R[3][1]=a-beta*V[0];	
	R[3][2]=-beta*V[1];	
	R[3][3]=-beta*V[2];
	R[4][0]=gamma_s+V[0]*a;
	R[4][1]=-a-beta*V[0];	
	R[4][2]=-beta*V[1];	
	R[4][3]=-beta*V[2];
	R[3][4]=R[4][4]=beta;	

	mtx_mult(L,R,Jacob);

	if (signal>0) {
		Q[0]=left.rho;
		Q[1]=left.rho*left.Vn[0];
		Q[2]=left.rho*left.Vn[1];
		Q[3]=left.rho*left.Vn[2];
		Q[4]=left.rho*left.H-left.p;
	} else {
		Q[0]=right.rho;
		Q[1]=right.rho*right.Vn[0];
		Q[2]=right.rho*right.Vn[1];
		Q[3]=right.rho*right.Vn[2];
		Q[4]=right.rho*right.H-right.p;
	}

	for (int i=0;i<5;++i) for (int j=0;j<5;++j) fluxNormal[i]+=Jacob[i][j]*Q[j];
	
	} // End signal loop

	return;
}


void mtx_out(double L[5][5],int dim){

      //int  dim=L[0].size();
      for(int row=0;row<dim;row++){
         for(int col=0;col<dim;col++) cout<<"L["<<row<<"]["<<col<<"]="<<L[row][col]<<'\t';
         cout<<endl;
      }
}
