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

void NavierStokes::mpi_init(void) {
	
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	return;
}

void NavierStokes::mpi_update_ghost_primitives(void) {
	
	for (int g=0;g<grid[gid].ghostCount;++g) {
		update[0].ghost(g)=p.ghost(g);
		update[1].ghost(g)=V.ghost(g)[0];
		update[2].ghost(g)=V.ghost(g)[1];
		update[3].ghost(g)=V.ghost(g)[2];
		update[4].ghost(g)=T.ghost(g);		
	}
	
	// The Following is convenient but not efficient
	/*
	p.mpi_update();
	V.mpi_update();
	T.mpi_update();
	 */
	int id;
	vector<double> sendBuffer;
	vector<double> recvBuffer;

	long unsigned int max_send_size=0;
	long unsigned int max_recv_size=0;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			max_send_size=max(max_send_size,5*grid[gid].sendCells[proc].size());
			max_recv_size=max(max_recv_size,5*grid[gid].recvCells[proc].size());
		}
	}
	sendBuffer.reserve(max_send_size); recvBuffer.reserve(max_recv_size);
			
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			sendBuffer.resize(5*grid[gid].sendCells[proc].size());
			recvBuffer.resize(5*grid[gid].recvCells[proc].size());
			
			for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
				id=grid[gid].maps.cellGlobal2Local[grid[gid].sendCells[proc][g]];
				sendBuffer[g*5]=p.cell(id);
				sendBuffer[g*5+1]=V.cell(id)[0];
				sendBuffer[g*5+2]=V.cell(id)[1];
				sendBuffer[g*5+3]=V.cell(id)[2];
				sendBuffer[g*5+4]=T.cell(id);
			}
			
			MPI_Sendrecv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,proc,0,&recvBuffer[0],recvBuffer.size(),MPI_DOUBLE,proc,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].maps.ghostGlobal2Local[grid[gid].recvCells[proc][g]];
				p.ghost(id)=recvBuffer[g*5];
				V.ghost(id)[0]=recvBuffer[g*5+1];
				V.ghost(id)[1]=recvBuffer[g*5+2];
				V.ghost(id)[2]=recvBuffer[g*5+3];
				T.ghost(id)=recvBuffer[g*5+4];				
			}
		}
	}
	
	sendBuffer.clear(); recvBuffer.clear();
	
	
	for (int g=0;g<grid[gid].ghostCount;++g) {
		update[0].ghost(g)=p.ghost(g)-update[0].ghost(g);
		update[1].ghost(g)=V.ghost(g)[0]-update[1].ghost(g);
		update[2].ghost(g)=V.ghost(g)[1]-update[2].ghost(g);
		update[3].ghost(g)=V.ghost(g)[2]-update[3].ghost(g);
		update[4].ghost(g)=T.ghost(g)-update[4].ghost(g);
		rho.ghost(g)=material.rho(p.ghost(g),T.ghost(g));
	}
	
	return;
} 

void NavierStokes::mpi_update_ghost_gradients(void) {
	
	// The Following is convenient but not efficient
	/*
	gradp.mpi_update();
	gradu.mpi_update();
	gradv.mpi_update();
	gradw.mpi_update();
	gradT.mpi_update();
	gradrho.mpi_update();
	*/
	
	int id;
	vector<double> sendBuffer;
	vector<double> recvBuffer;
	
	long unsigned int max_send_size=0;
	long unsigned int max_recv_size=0;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			max_send_size=max(max_send_size,15*grid[gid].sendCells[proc].size());
			max_recv_size=max(max_recv_size,15*grid[gid].recvCells[proc].size());
		}
	}
	sendBuffer.reserve(max_send_size); recvBuffer.reserve(max_recv_size);
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			sendBuffer.resize(15*grid[gid].sendCells[proc].size());
			recvBuffer.resize(15*grid[gid].recvCells[proc].size());
			
			for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
				id=grid[gid].maps.cellGlobal2Local[grid[gid].sendCells[proc][g]];
				for (int i=0;i<3;++i) {
					sendBuffer[g*15+i]=gradp.cell(id)[i];
					sendBuffer[g*15+3+i]=gradu.cell(id)[i];
					sendBuffer[g*15+6+i]=gradv.cell(id)[i];
					sendBuffer[g*15+9+i]=gradw.cell(id)[i];
					sendBuffer[g*15+12+i]=gradT.cell(id)[i];
				}
			}
			
			MPI_Sendrecv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,proc,0,&recvBuffer[0],recvBuffer.size(),MPI_DOUBLE,proc,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].maps.ghostGlobal2Local[grid[gid].recvCells[proc][g]];
				for (int i=0;i<3;++i) {
					gradp.ghost(id)[i]=recvBuffer[g*15+i];
					gradu.ghost(id)[i]=recvBuffer[g*15+3+i];
					gradv.ghost(id)[i]=recvBuffer[g*15+6+i];
					gradw.ghost(id)[i]=recvBuffer[g*15+9+i];
					gradT.ghost(id)[i]=recvBuffer[g*15+12+i];
				}				
			}
		}
	}
	
	sendBuffer.clear(); recvBuffer.clear();
	
	return;
} 

