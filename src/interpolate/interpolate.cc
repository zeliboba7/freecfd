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
#include "interpolate.h"

void Interpolate::init(void) {
	
	a3.resize(3); b3.resize(3); weights3.resize(3);
	for (int i=0;i<3;++i) a3[i].resize(3);
	
	a4.resize(4); b4.resize(4); weights4.resize(4);
	for (int i=0;i<4;++i) a4[i].resize(4);
	
	return;
}

void Interpolate::calculate_weights(bool is_internal) {
	
	// Sort the stencil based on closest distance
	// and truncate based on max_stencil_size
	sort_stencil(is_internal);
	
	weights.resize(stencil.size());
	for (int i=0;i<weights.size();++i) weights[i]=0.;
		 
	if (interpolate_tetra())     kind=4;
	else if (interpolate_tri())  kind=3;
	else if (interpolate_line()) kind=2;
	else { interpolate_point();  kind=1; }
	
	return;
}

void Interpolate::flush(void) {
	
	stencil.clear();
	stencil_indices.clear();
	weights.clear();
	
	return;
}

void Interpolate::sort_stencil(bool is_internal) {

	int stencil_size=stencil.size();
	vector<double> distances;
	vector<int> order (min(stencil_size,max_stencil_size),0);
	vector<double> min_distances (order.size(),1.e20);
	
	// Store distances to the point
	for (vit=stencil.begin();vit!=stencil.end();vit++) distances.push_back(fabs(*vit-point));
	
	vector<Vec3D> new_stencil (order.size(),0.);
	vector<int> new_stencil_indices (order.size(),0);

	//is_internal=false;
	
	// Sort distances
	if (is_internal) { // Keep the first 2 (parent and neighbor) in the stencil unchanged
		order[0]=0;
		order[1]=1;
		for (int d=2;d<distances.size();++d) {
			for (int m=2;m<min_distances.size();++m) {
				if (distances[d]<=min_distances[m]) {
					for (int i=order.size()-1;i>m;--i) order[i]=order[i-1];
					for (int i=order.size()-1;i>m;--i) min_distances[i]=min_distances[i-1];
					order[m]=d;
					min_distances[m]=distances[d];
					break;
				}
			}
		}
	} else {
		for (int d=0;d<distances.size();++d) {
			for (int m=0;m<min_distances.size();++m) {
				if (distances[d]<=min_distances[m]) {
					for (int i=order.size()-1;i>m;--i) order[i]=order[i-1];
					for (int i=order.size()-1;i>m;--i) min_distances[i]=min_distances[i-1];
					order[m]=d;
					min_distances[m]=distances[d];
					break;
				}
			}
		}

	}

	for (int i=0;i<order.size();++i) {
		new_stencil[i]=stencil[order[i]];
		new_stencil_indices[i]=stencil_indices[order[i]];
	}

	stencil=new_stencil;
	stencil_indices=new_stencil_indices;
	
	distances.clear();
	order.clear();
	min_distances.clear();
	new_stencil.clear();
	
	return;
}

bool Interpolate::interpolate_tetra(void) {
	
	if (dimension!=3) return false;
	if (stencil.size()<4) return false;
	
	for (int i=0;i<weights.size();++i) weights[i]=0.;
	
	bool is_success=false;
	
	for (int s1=0;s1<stencil.size()-2;++s1) {
		for (int s2=s1+1;s2<stencil.size()-1;++s2) {
			for (int s3=s2+1;s3<stencil.size();++s3) {
				for (int s4=s3+1;s4<stencil.size();++s4) {
					// Evaluate viability of the constructed tetra for interpolation
					edge1=stencil[s2]-stencil[s1];
					edge2=stencil[s3]-stencil[s1];
					edge3=stencil[s3]-stencil[s2];
					edge4=stencil[s4]-stencil[s1];
					edge5=stencil[s4]-stencil[s2];
					edge6=stencil[s4]-stencil[s3];
					ave_edge=(fabs(edge1)+fabs(edge2)+fabs(edge3)+fabs(edge4)+fabs(edge5)+fabs(edge6))/6.;
					// Area and normal of the base triangle (formed by edges 1,2,3)
					planeNormal=edge1.cross(edge2);
					area=0.5*fabs(planeNormal);
					planeNormal=planeNormal.norm();
					volume=area*fabs(edge4.dot(planeNormal))/3.;
					equilateral_volume=sqrt(2.)/12.*ave_edge*ave_edge*ave_edge;		
					skewness=1.-volume/equilateral_volume; // really skewed -> 1, perfect -> 0
					// If the tetra is too skewed, skip it
					if (skewness>(skewness_tolerance*skewness_tolerance*skewness_tolerance)) continue;
					// At this point, there is a decent tetra
					is_success=true;
					// Weight its contribution by how close its centroid is to the point to interpolate
					centroid=(stencil[s1]+stencil[s2]+stencil[s3]+stencil[s4])/4.;
					distance=fabs(point-centroid)/ave_edge; // normalized distance
					distance=max(distance,1.e-2); // limit minimum in case point and centroid coincide
					tetra_weight=1./(distance*distance); // far away -> 0 , really close >> 1
					// Now to the actual interpolation
					basis1=edge1;
					basis1=basis1.norm();
					basis2=-1.*(basis1.cross(planeNormal)).norm();
					basis3=(basis1.cross(basis2)).norm();

					p_point-=stencil[s1]; // stencil[s1] is acting as the origin of a new coordinate system
					
					// Form the linear system
					a4[0][0]=0.;
					a4[0][1]=edge1.dot(basis1);
					a4[0][2]=edge2.dot(basis1);
					a4[0][3]=edge4.dot(basis1);
					
					a4[1][0]=0.;
					a4[1][1]=edge1.dot(basis2);
					a4[1][2]=edge2.dot(basis2);
					a4[1][2]=edge4.dot(basis2);
					
					a4[2][0]=0.;
					a4[2][1]=edge1.dot(basis3);
					a4[2][2]=edge2.dot(basis3);
					a4[2][2]=edge4.dot(basis3);
					
					a4[3][0]=1.;
					a4[3][1]=1.;
					a4[3][2]=1.;
					a4[3][2]=1.;
					
					b4[0]=p_point.dot(basis1);
					b4[1]=p_point.dot(basis2);
					b4[2]=p_point.dot(basis3);
					b4[3]=1.;
					
					// Solve the 4x4 linear system by Gaussion Elimination
					gelimd(a4,b4,weights4);
					// Let's see if the linear system solution is good
					weightSum=0.;
					for (int i=0;i<4;++i) weightSum+=weights4[i];
					if (fabs(weightSum-1.)>1.e-8) {
					//	cout << "[W rank=" << Rank << "] Tetra interpolation weightSum=" << setprecision(8) << weightSum << " is not unity" << endl;
					//	cout << "[W rank=" << Rank << "] Switching to tri interpolation" << endl;
						return false;
					}
					
					weights[s1]+=tetra_weight*weights4[0];
					weights[s2]+=tetra_weight*weights4[1];
					weights[s3]+=tetra_weight*weights4[2];
					weights[s4]+=tetra_weight*weights4[3];
					
				}
			}
		}
	}
	
	weightSum=0.;
	for (int i=0;i<weights.size();++i) weightSum+=weights[i];
	
	for (int i=0;i<weights.size();++i) weights[i]/=weightSum;
	
	return is_success;
	
}

bool Interpolate::interpolate_tri(void) {
	
	if (dimension<2) return false;
	if (stencil.size()<3) return false;
	
	for (int i=0;i<weights.size();++i) weights[i]=0.;
	
	bool is_success=false;

	for (int s1=0;s1<stencil.size()-2;++s1) {
		for (int s2=s1+1;s2<stencil.size()-1;++s2) {
			for (int s3=s2+1;s3<stencil.size();++s3) {
				// Evaluate viability of the constructed triangle for interpolation
				edge1=stencil[s2]-stencil[s1];
				edge2=stencil[s3]-stencil[s1];
				edge3=stencil[s3]-stencil[s2];
				planeNormal=edge1.cross(edge2);
				area=0.5*fabs(planeNormal);
				planeNormal=planeNormal.norm();
				ave_edge=(fabs(edge1)+fabs(edge2)+fabs(edge3))/3.;
				equilateral_area=0.433*ave_edge*ave_edge;
				skewness=1.-area/equilateral_area; // really skewed -> 1, perfect -> 0
				// If the triangle is too skewed, skip it
				if (skewness>skewness_tolerance*skewness_tolerance) continue;
				// At this point, there is a decent triangle
				is_success=true;
				// Weight its contribution by how close its centroid is to the point to interpolate
				// First project the point to the plane of the triangle
				p_point=point-point.dot(planeNormal)*planeNormal;
				centroid=(stencil[s1]+stencil[s2]+stencil[s3])/3.;
				distance=fabs(p_point-centroid)/ave_edge; // normalized distance
				distance=max(distance,1.e-2); // limit minimum in case point and centroid coincide
				tri_weight=1./(distance*distance); // far away -> 0 , really close >> 1
				// Now to the actual interpolation
				basis1=edge1;
				basis1=basis1.norm();
				basis2=-1.*(basis1.cross(planeNormal)).norm();
				p_point-=stencil[s1]; // stencil[s1] is acting as the origin of a new coordinate system
				
				// Form the linear system
				a3[0][0]=0.;
				a3[0][1]=edge1.dot(basis1);
				a3[0][2]=edge2.dot(basis1);
				a3[1][0]=0.;
				a3[1][1]=edge1.dot(basis2);
				a3[1][2]=edge2.dot(basis2);
				a3[2][0]=1.;
				a3[2][1]=1.;
				a3[2][2]=1.;
				b3[0]=p_point.dot(basis1);
				b3[1]=p_point.dot(basis2);
				b3[2]=1.;
				
				// Solve the 3x3 linear system by Gaussion Elimination
				gelimd(a3,b3,weights3);
				// Let's see if the linear system solution is good
				weightSum=0.;
				for (int i=0;i<3;++i) weightSum+=weights3[i];
				if (fabs(weightSum-1.)>1.e-8) {
					cout << "[W rank=" << Rank << "] Tri interpolation weightSum=" << setprecision(8) << weightSum << " is not unity" << endl;
					cout << "[W rank=" << Rank << "] Switching to line interpolation" << endl;
					return false;
				}
				
				weights[s1]+=tri_weight*weights3[0];
				weights[s2]+=tri_weight*weights3[1];
				weights[s3]+=tri_weight*weights3[2];

			}
		}
	}
	
	weightSum=0.;
	for (int i=0;i<weights.size();++i) weightSum+=weights[i];
	
	for (int i=0;i<weights.size();++i) weights[i]/=weightSum;

	return is_success;

}

bool Interpolate::interpolate_line(void) {
	
	if (stencil.size()<2) return false;
	
	edge1=stencil[1]-stencil[0]; // line vector
	basis1=edge1; 
	basis1=basis1.norm(); // direction of the line
	
	// Try this ! than fix above if it works properly
	//cout << edge1 << "\t" << edge1.norm() << "\t" << edge1 << endl;
	
	// point moved to new coordinate system
	p_point=point-stencil[0]; 
	// point projected on the line
	p_point=p_point.dot(basis1)*basis1;

	weights[1]=p_point.dot(basis1)/fabs(edge1);
	weights[0]=1.-weights[1];
	
	weights.resize(2);
	stencil_indices.resize(2);

	return true;
}

bool Interpolate::interpolate_point(void) {
	
	weights.resize(1); 
	stencil_indices.resize(1);
	weights[0]=1.;
	
	return true;
}
