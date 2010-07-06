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
#include "grid.h"
#include "inputs.h"
#include "bc.h"
#include "bc_interface.h"
#include "ns.h"
#include "hc.h"
#include "commons.h"

extern InputFile input;
extern vector<Grid> grid;
extern vector<vector<BCregion> > bc;
extern vector<vector<BC_Interface> > interface;

extern vector<NavierStokes> ns;
extern vector<HeatConduction> hc;


void bc_interface_sync(void) {
	
	// Loop all the interfaces 
	for (int gid=0;gid<grid.size();++gid) {
		for (int i=0;i<interface[gid].size();++i) {
			int Rank=grid[gid].Rank;
			int np=grid[gid].np;
			int donor_bc=interface[gid][i].donor_bc;
			int donor_grid=interface[gid][i].donor_grid;
			int	recvSize=grid[donor_grid].globalBoundaryFaceCount[donor_bc];
			vector<double> sendBuffer;
			vector<double> recvBuffer (recvSize,0.);
			vector<int> recvSizes (np,0);
			vector<int> displs (np,0);
			for (int p=0;p<np;++p) {
				recvSizes[p]=grid[donor_grid].boundaryFaceCount[donor_bc][p];
				if (p>0) displs[p]=recvSizes[p-1];
			}
			int counter=0;
			for (int f=0;f<grid[donor_grid].faceCount;++f) {
				if (grid[donor_grid].face[f].bc==donor_bc) {
					//cout << ns[donor_grid].qdot.face(f) << endl;
					// TODO: generalize by adding and equation int to the bc_interface class
					//cout << interface[gid][i].donor_var << "\t" << ns[donor_grid].qdot.face(f) << endl;
					if (interface[gid][i].donor_var=="T") {
						if (interface[gid][i].donor_eqn==NS) sendBuffer.push_back(ns[donor_grid].T.face(f));
						else if (interface[gid][i].donor_eqn==HEAT) sendBuffer.push_back(hc[donor_grid].T.face(f));
					}
					else if (interface[gid][i].donor_var=="qdot") {
						//if (interface[gid][i].donor_eqn==NS) sendBuffer.push_back(ns[donor_grid].qdot.bc(interface[gid][i].donor_bc,f));
						if (interface[gid][i].donor_eqn==NS) sendBuffer.push_back(ns[donor_grid].qdot.face(f));
						else if (interface[gid][i].donor_eqn==HEAT) sendBuffer.push_back(hc[donor_grid].qdot.face(f));
					}
					counter++;
				}
			}

			MPI_Allgatherv(&sendBuffer[0],sendBuffer.size(),MPI_DOUBLE,&recvBuffer[0],&recvSizes[0],&displs[0],MPI_DOUBLE,MPI_COMM_WORLD);

			for (int j=0;j<interface[gid][i].donor_point.size();++j) interface[gid][i].donor_data[j]=recvBuffer[j];
			
		}

	}

	for (int gid=0;gid<grid.size();++gid) {
		int count=0;
		for (int i=0;i<interface[gid].size();++i) {
			for (int f=0;f<grid[gid].faceCount;++f) {
				if (grid[gid].face[f].bc==interface[gid][i].recv_bc) {
					if (interface[gid][i].donor_var=="T") {
						if (interface[gid][i].recv_eqn==NS) ns[gid].T.face(f)=interface[gid][i].donor_data[interface[gid][i].donor_index[count]];
						else if (interface[gid][i].recv_eqn==HEAT) hc[gid].T.face(f)=interface[gid][i].donor_data[interface[gid][i].donor_index[count]];
					}
					else if (interface[gid][i].donor_var=="qdot") {
						if (interface[gid][i].recv_eqn==NS) ns[gid].qdot.face(f)=interface[gid][i].donor_data[interface[gid][i].donor_index[count]];
						else if (interface[gid][i].recv_eqn==HEAT) hc[gid].qdot.face(f)=interface[gid][i].donor_data[interface[gid][i].donor_index[count]];
					}
					count++;
				}
			}
		}
	}
	
	
	return;
}
