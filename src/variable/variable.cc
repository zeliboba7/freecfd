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

#include "variable.h"

template <>
void Variable<double>::mpi_update (void) {
	
	int np=grid[gid].np;
	int Rank=grid[gid].Rank;
	int id;
	
	vector<double> sendBuffer;
	vector<double> recvBuffer;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			sendBuffer.resize(grid[gid].sendCells[proc].size());
			recvBuffer.resize(grid[gid].recvCells[proc].size());
			
			for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
				id=grid[gid].maps.cellGlobal2Local[grid[gid].sendCells[proc][g]];
				sendBuffer[g]=cell(id);
			}
			
			MPI_Sendrecv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,proc,0,&recvBuffer[0],recvBuffer.size(),MPI_DOUBLE,proc,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].maps.ghostGlobal2Local[grid[gid].recvCells[proc][g]];
				ghost(id)=recvBuffer[g];
			}
		}
	}

	sendBuffer.clear(); recvBuffer.clear();
	
	return;
}

template <>
void Variable<Vec3D>::mpi_update (void) {
	
	int np=grid[gid].np;
	int Rank=grid[gid].Rank;
	int id;
	
	vector<double> sendBuffer;
	vector<double> recvBuffer;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			sendBuffer.resize(3*grid[gid].sendCells[proc].size());
			recvBuffer.resize(3*grid[gid].recvCells[proc].size());
		
			for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
				id=grid[gid].maps.cellGlobal2Local[grid[gid].sendCells[proc][g]];
				for (int i=0;i<3;++i) sendBuffer[g*3+i]=cell(id)[i];
			}
			
			MPI_Sendrecv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,proc,0,&recvBuffer[0],recvBuffer.size(),MPI_DOUBLE,proc,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].maps.ghostGlobal2Local[grid[gid].recvCells[proc][g]];
				for (int i=0;i<3;++i) ghost(id)[i]=recvBuffer[g*3+i];
			}
		}
	
	}

	sendBuffer.clear(); recvBuffer.clear();

	return;
}

template <>
double Variable<double>::cell2node (int c, int n) {
	Vec3D grad;
	grad=cell_gradient(c);
	return cell(c)+grad.dot(grid[gid].node[n]-grid[gid].cell[c].centroid);
}

template <>
Vec3D Variable<Vec3D>::cell2node (int c, int n) {
	vector<Vec3D> grad (3,0.);
	Vec3D result;
	grad=cell_gradient(c);
	result[0]=grad[0].dot(grid[gid].node[n]-grid[gid].cell[c].centroid);
	result[1]=grad[1].dot(grid[gid].node[n]-grid[gid].cell[c].centroid);
	result[2]=grad[2].dot(grid[gid].node[n]-grid[gid].cell[c].centroid);
	return cell(c)+result;
}







