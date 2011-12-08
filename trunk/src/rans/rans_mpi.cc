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
		
	long unsigned int max_send_size=0;
	long unsigned int max_recv_size=0;
	
	mpi_send_offset.resize(np);
	mpi_recv_offset.resize(np);
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			mpi_send_offset[proc]=max_send_size;
			mpi_recv_offset[proc]=max_recv_size;
			max_send_size+=6*grid[gid].sendCells[proc].size();
			max_recv_size+=6*grid[gid].recvCells[proc].size();
		}
	}
	sendBuffer.resize(max_send_size); recvBuffer.resize(max_recv_size);
	
	send_req_count=0;
	recv_req_count=0;
	
	for (int proc=0;proc<np;++proc) { 
		if (Rank!=proc) {
			if (grid[gid].sendCells[proc].size()!=0) send_req_count++;;
			if (grid[gid].recvCells[proc].size()!=0) recv_req_count++;
		}
	}
	
	send_request.resize(send_req_count);
	recv_request.resize(recv_req_count);
	
	return;
}

void RANS::mpi_update_ghost_primitives(void) {

	for (int g=grid[gid].partition_ghosts_begin;g<=grid[gid].partition_ghosts_end;++g) {
		update[0].cell(g)=k.cell(g);
		update[1].cell(g)=omega.cell(g);
	}
	
	int id,offset;
	
	send_req_count=0; recv_req_count=0;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			if (grid[gid].sendCells[proc].size()!=0) {
				for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
					id=grid[gid].sendCells[proc][g];
					offset=mpi_send_offset[proc];
					sendBuffer[offset+g*2]=k.cell(id);
					sendBuffer[offset+g*2+1]=omega.cell(id);
				}
				
				MPI_Isend(&sendBuffer[offset],grid[gid].sendCells[proc].size()*2,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&send_request[send_req_count]);
				send_req_count++;
			}
			
			if (grid[gid].recvCells[proc].size()!=0) {
				offset=mpi_recv_offset[proc];
				MPI_Irecv(&recvBuffer[offset],grid[gid].recvCells[proc].size()*2,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&recv_request[recv_req_count]);
				recv_req_count++;
			}
		}
	}
	
    MPI_Waitall(recv_request.size(),&recv_request[0],MPI_STATUS_IGNORE);
	
	for (int proc=0;proc<np;++proc) { 
		if (Rank!=proc) {
			offset=mpi_recv_offset[proc];
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].recvCells[proc][g];
				k.cell(id)=recvBuffer[offset+g*2];
				omega.cell(id)=recvBuffer[offset+g*2+1];
			}
		}
	}

	for (int g=grid[gid].partition_ghosts_begin;g<=grid[gid].partition_ghosts_end;++g) {
		update[0].cell(g)=k.cell(g)-update[0].cell(g);
		update[1].cell(g)=omega.cell(g)-update[1].cell(g);
	}
	
	return;
} 

void RANS::mpi_update_ghost_gradients(void) {

	int id,offset;
	
	send_req_count=0; recv_req_count=0;
	
	for (int proc=0;proc<np;++proc) {
		if (Rank!=proc) {
			if (grid[gid].sendCells[proc].size()!=0) {
				for (int g=0;g<grid[gid].sendCells[proc].size();++g) {
					id=grid[gid].sendCells[proc][g];
					offset=mpi_send_offset[proc];
					for (int i=0;i<3;++i) {
						sendBuffer[offset+g*6+i]=gradk.cell(id)[i];
						sendBuffer[offset+g*6+3+i]=gradomega.cell(id)[i];
					}
				}
				
				MPI_Isend(&sendBuffer[offset],grid[gid].sendCells[proc].size()*6,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&send_request[send_req_count]);
				send_req_count++;
			}
			
			if (grid[gid].recvCells[proc].size()!=0) {
				offset=mpi_recv_offset[proc];
				MPI_Irecv(&recvBuffer[offset],grid[gid].recvCells[proc].size()*6,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&recv_request[recv_req_count]);
				recv_req_count++;
			}
		}
	}
	
	MPI_Waitall(recv_request.size(),&recv_request[0],MPI_STATUS_IGNORE);
	
	for (int proc=0;proc<np;++proc) { 
		if (Rank!=proc) {
			offset=mpi_recv_offset[proc];
			for (int g=0;g<grid[gid].recvCells[proc].size();++g) {
				id=grid[gid].recvCells[proc][g];
				for (int i=0;i<3;++i) {
					gradk.cell(id)[i]=recvBuffer[offset+g*6+i];
					gradomega.cell(id)[i]=recvBuffer[offset+g*6+3+i];
				}
			}
		}
	}

	return;
} 

