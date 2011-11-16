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

extern void roe_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Gamma,double &weightL);
extern void vanLeer_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Gamma,double Pref,double &weightL);
extern void AUSMplusUP_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Gamma,double Pref,double Minf,double &weightL);
extern void SD_SLAU_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Pref,double &weightL);
extern void Stegger_Warming_flux(NS_Cell_State &left,NS_Cell_State &right,double diss_factor,double closest_wall_distance,double wdiss,double bl_height,MATERIAL &material,double fluxNormal[],double &weightL);
void flux_from_right(NS_Cell_State &right,double fluxNormal[]);

void NavierStokes::convective_face_flux(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,double flux[]) {

	double fluxNormal[5];

	if (face.bc>=0 && bc[gid][face.bc].type==INLET) {
		flux_from_right(right,fluxNormal);
		weightL.face(face.index)=0.;
	} else if (convective_flux_function==ROE) {
		roe_flux(left,right,fluxNormal,material.gamma,weightL.face(face.index));
	} else if (convective_flux_function==VAN_LEER) {
		vanLeer_flux(left,right,fluxNormal,material.gamma,material.Pref,weightL.face(face.index));
	} else if (convective_flux_function==AUSM_PLUS_UP) {
		AUSMplusUP_flux(left,right,fluxNormal,material.gamma,material.Pref,Minf,weightL.face(face.index));
	} else if (convective_flux_function==SD_SLAU) {
		SD_SLAU_flux(left,right,fluxNormal,material.Pref,weightL.face(face.index));
	} else if (convective_flux_function==SW) {
		Stegger_Warming_flux(left,right,
					grid[gid].face[face.index].dissipation_factor,
					grid[gid].face[face.index].closest_wall_distance,
					wdiss,bl_height,
					material,fluxNormal,weightL.face(face.index));
	} 
	
	flux[0] = fluxNormal[0]*face.area;
	flux[1] = (fluxNormal[1]*face.normal[0]+fluxNormal[2]*face.tangent1[0]+fluxNormal[3]*face.tangent2[0])*face.area;
	flux[2] = (fluxNormal[1]*face.normal[1]+fluxNormal[2]*face.tangent1[1]+fluxNormal[3]*face.tangent2[1])*face.area;
	flux[3] = (fluxNormal[1]*face.normal[2]+fluxNormal[2]*face.tangent1[2]+fluxNormal[3]*face.tangent2[2])*face.area;
	flux[4] = fluxNormal[4]*face.area;
	mdot.face(face.index)=fluxNormal[0];

	return;
} // end face flux

void flux_from_right(NS_Cell_State &right,double fluxNormal[]) {

	double mdot=right.rho*right.Vn[0];

	fluxNormal[0]=mdot;
	fluxNormal[1]=mdot*right.Vn[0]+right.p;
	fluxNormal[2]=mdot*right.Vn[1];
	fluxNormal[3]=mdot*right.Vn[2];
	fluxNormal[4]=mdot*right.H;

	return;
}
