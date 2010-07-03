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
#include "bc_interface.h"

void BC_Interface::setup(void) {
	
	int Rank=grid[recv_grid].Rank;
	int np=grid[recv_grid].np;
	
	int	recvSize=3*grid[donor_grid].globalBoundaryFaceCount[donor_bc];
	donor_data.resize(recvSize/3);
	donor_point.resize(recvSize/3);
	vector<double> sendBuffer;
	vector<double> recvBuffer (recvSize,0.);
	
	vector<int> recvSizes (np,0);
	vector<int> displs (np,0);
	for (int p=0;p<np;++p) {
		recvSizes[p]=grid[donor_grid].boundaryFaceCount[donor_bc][p];
		if (p>0) displs[p]=recvSizes[p-1];
	}

	for (int f=0;f<grid[donor_grid].faceCount;++f) {
		if (grid[donor_grid].face[f].bc==donor_bc) {
			sendBuffer.push_back(grid[donor_grid].face[f].centroid[0]);
			sendBuffer.push_back(grid[donor_grid].face[f].centroid[1]);
			sendBuffer.push_back(grid[donor_grid].face[f].centroid[2]);
		}
	}

	MPI_Allgatherv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,&recvBuffer[0],&recvSizes[0],&displs[0],MPI_DOUBLE,MPI_COMM_WORLD);

	
	for (int i=0;i<donor_point.size();++i) {
		donor_point[i][0]=recvBuffer[3*i];
		donor_point[i][1]=recvBuffer[3*i+1];
		donor_point[i][2]=recvBuffer[3*i+2];
	}
	
	// TODO: For now, establish only single closest point search
	// In the future, add LTI
	donor_index.resize(grid[recv_grid].boundaryFaceCount[recv_bc][Rank]);
	double min_distance,distance;
	int i_closest;
	int counter=0;
	for (int f=0;f<grid[recv_grid].faceCount;++f) {
		if (grid[recv_grid].face[f].bc==recv_bc) {
			// Search for the closest point in the donor point cloud
			min_distance=1.e20;
			for (int i=0;i<donor_point.size();++i) {
				distance=fabs(donor_point[i]-grid[recv_grid].face[f].centroid);
				if (distance<min_distance) {
					min_distance=distance;
					i_closest=i;
				}
			}
			donor_index[counter]=i_closest;
			counter++;
		}
	}

	return;
}
