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

void RANS::mpi_init(void) {
	
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	return;
}

void RANS::mpi_update_ghost_primitives(void) {
	
	for (int g=0;g<grid[gid].ghostCount;++g) {
		update[0].ghost(g)=k.ghost(g);
		update[1].ghost(g)=omega.ghost(g);
	}
	
	/*
	k.mpi_update();
	omega.mpi_update();
	mu_t.mpi_update();
	*/
	
	int id;
	vector<double> sendBuffer;
	vector<double> recvBuffer;
	
	long unsigned int max_send_size=0;
	long unsigned int max_recv_size=0;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			max_send_size=max(max_send_size,3*grid[gid].sendCells[proc].size());
			max_recv_size=max(max_recv_size,3*grid[gid].recvCells[proc].size());
		}
	}
	sendBuffer.reserve(max_send_size); recvBuffer.reserve(max_recv_size);
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			sendBuffer.resize(3*grid[gid].sendCells[proc].size());
			recvBuffer.resize(3*grid[gid].recvCells[proc].size());
			
			for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
				id=grid[gid].maps.cellGlobal2Local[grid[gid].sendCells[proc][g]];
				sendBuffer[g*3]=k.cell(id);
				sendBuffer[g*3+1]=omega.cell(id);
				sendBuffer[g*3+2]=mu_t.cell(id);
			}
			
			MPI_Sendrecv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,proc,0,&recvBuffer[0],recvBuffer.size(),MPI_DOUBLE,proc,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].maps.ghostGlobal2Local[grid[gid].recvCells[proc][g]];
				k.ghost(id)=recvBuffer[g*3];
				omega.ghost(id)=recvBuffer[g*3+1];
				mu_t.ghost(id)=recvBuffer[g*3+2];			
			}
		}
	}
	

	for (int g=0;g<grid[gid].ghostCount;++g) {
		update[0].ghost(g)=k.ghost(g)-update[0].ghost(g);
		update[1].ghost(g)=omega.ghost(g)-update[1].ghost(g);
	}
	
	return;
} 

void RANS::mpi_update_ghost_gradients(void) {
	
	/*
	gradk.mpi_update();
	gradomega.mpi_update();
	*/
	
	int id;
	vector<double> sendBuffer;
	vector<double> recvBuffer;
	
	long unsigned int max_send_size=0;
	long unsigned int max_recv_size=0;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			max_send_size=max(max_send_size,6*grid[gid].sendCells[proc].size());
			max_recv_size=max(max_recv_size,6*grid[gid].recvCells[proc].size());
		}
	}
	sendBuffer.reserve(max_send_size); recvBuffer.reserve(max_recv_size);
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			sendBuffer.resize(6*grid[gid].sendCells[proc].size());
			recvBuffer.resize(6*grid[gid].recvCells[proc].size());
			
			for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
				id=grid[gid].maps.cellGlobal2Local[grid[gid].sendCells[proc][g]];
				for (int i=0;i<3;++i) {
					sendBuffer[g*6+i]=gradk.cell(id)[i];
					sendBuffer[g*6+3+i]=gradomega.cell(id)[i];
				}
			}
			
			MPI_Sendrecv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,proc,0,&recvBuffer[0],recvBuffer.size(),MPI_DOUBLE,proc,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].maps.ghostGlobal2Local[grid[gid].recvCells[proc][g]];
				for (int i=0;i<3;++i) {
					gradk.ghost(id)[i]=recvBuffer[g*6+i];
					gradomega.ghost(id)[i]=recvBuffer[g*6+3+i];
				}				
			}
		}
	}
	
	sendBuffer.clear(); recvBuffer.clear();

	return;
} 

