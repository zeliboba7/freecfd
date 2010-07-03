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

string int2str(int number) ;

Grid::Grid() {
	// Just get the current processor's rank
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	// And total number of processors
	MPI_Comm_size(MPI_COMM_WORLD, &np);
}

void Grid::read(string fname) {
	fstream file;
	fileName=fname;
	file.open(fileName.c_str());
	if (file.is_open()) {
		if (Rank==0) cout << "[I] Found grid file " << fileName  << endl;
		file.close();		
		readCGNS();
	} else {
		if (Rank==0) cerr << "[E] Grid file "<< fileName << " could not be found." << endl;
		exit(1);
	}
	return;
}

void Grid::setup(void) {
	partition();
	create_nodes_cells();
	mesh2dual();
	create_faces();
	create_ghosts();
	trim_memory();
	areas_volumes();
	mpi_handshake();
	mpi_get_ghost_geometry();
	return;
}
	
void Grid::trim_memory() {
	// A trick for shrinking vector capacities to just the right sizes
	
	for (int n=0;n<nodeCount;++n) {
		vector<int> (node[n].cells).swap(node[n].cells);
		vector<int> (node[n].ghosts).swap(node[n].ghosts);
		vector<int> (node[n].faces).swap(node[n].faces);
	}
	
	for (int f=0;f<faceCount;++f) { 
		vector<int> (face[f].nodes).swap(face[f].nodes);
	}

	for (int c=0;c<cellCount;++c) { 
		vector<int> (cell[c].nodes).swap(cell[c].nodes);
		vector<int> (cell[c].faces).swap(cell[c].faces);
		vector<int> (cell[c].neighborCells).swap(cell[c].neighborCells);
		vector<int> (cell[c].ghosts).swap(cell[c].ghosts);
	}
	
	for (int g=0;g<ghostCount;++g) { 
		vector<int> (ghost[g].cells).swap(ghost[g].cells);
	}
	
	vector<Node> (node).swap(node);
	vector<Face> (face).swap(face);
	vector<Cell> (cell).swap(cell);
	vector<Ghost> (ghost).swap(ghost);
	
	vector<int> (partitionOffset).swap(partitionOffset);
		
	// Destroy grid raw data
	raw.node.clear();
	raw.cellConnIndex.clear();
	raw.cellConnectivity.clear();
	// First save the number of boundary conditions
	bcCount=raw.bocoNodes.size();
	raw.bocoNodes.clear();
	raw.bocoNameMap.clear();
	
	return;
}

int Grid::translate(Vec3D begin, Vec3D end) {
	Vec3D diff;
	diff=end-begin;
	for (int n=0;n<globalNodeCount;++n) {
		raw.node[n]+=diff;
	}
	if (Rank==0) cout << "[I] Grid " << gid+1 << " is translated from " << begin << " to " << end << endl;
	return 1;
}

int Grid::scale(Vec3D anchor, double scale) {
	for (int n=0;n<globalNodeCount;++n) {
		raw.node[n]=anchor+scale*(raw.node[n]-anchor);
	}
	if (Rank==0) cout << "[I] Grid " << gid+1 << " is scaled by " << scale << " with anchor = " << anchor << endl;
	return 1;
}

int Grid::rotate(Vec3D anchor, Vec3D axis, double angle) {
	// Convert angle to radian
	angle*=4.*atan(1.)/180.;
	// Normalize axis;
	axis=axis.norm();
	Vec3D p;
	for (int n=0;n<globalNodeCount;++n) {
		p=raw.node[n];
		p-=anchor;
		raw.node[n][0]=axis[0]*(axis.dot(p))+(p[0]*(1.-axis[0]*axis[0])-axis[0]*(axis[1]*p[1]+axis[2]*p[2]))*cos(angle)+(-axis[2]*p[1]+axis[1]*p[2])*sin(angle);
		raw.node[n][1]=axis[1]*(axis.dot(p))+(p[1]*(1.-axis[1]*axis[1])-axis[1]*(axis[0]*p[0]+axis[2]*p[2]))*cos(angle)+(axis[2]*p[0]-axis[0]*p[2])*sin(angle);
		raw.node[n][2]=axis[2]*(axis.dot(p))+(p[2]*(1.-axis[2]*axis[2])-axis[2]*(axis[0]*p[0]+axis[1]*p[1]))*cos(angle)+(-axis[1]*p[0]+axis[0]*p[1])*sin(angle);
		raw.node[n]+=anchor;
	}
	angle*=180./(4.*atan(1.));
	if (Rank==0) cout << "[I] Grid " << gid+1 << " is rotated by " << angle << " degrees around axis = " << axis << " with anchor = " << anchor << endl;
	return 1;
}

int Grid::areas_volumes() {
	// NOTE The methodology here is generic hence slower than it could be:
	// i.e Even though volume/centroid of tetrahedra is easy enough, we still break it down to
	// smaller tetrahedras

	// Now loop through faces and calculate centroids and areas
	for (int f=0;f<faceCount;++f) {
		Vec3D centroid=0.;
		Vec3D areaVec=0.;
		Vec3D patchCentroid,patchArea;
		// Find an approxiamate centroid (true centroid for triangle)
		for (int n=0;n<face[f].nodeCount;++n) {
			centroid+=faceNode(f,n);
		}
		centroid/=double(face[f].nodeCount);
		
		// Sum the area as a patch of triangles formed by connecting two nodes and an interior point
		face[f].centroid=0.;
		areaVec=0.;
		int next;
		// First calculate face normal
		face[f].normal=(faceNode(f,2)-faceNode(f,1)).cross(faceNode(f,0)-faceNode(f,1));
		double length=fabs(face[f].normal);
		if (length<1.e-12) length=1.e-12;
		face[f].normal/=length;
		for (int n=0;n<face[f].nodeCount;++n) {
			next=n+1;
			if (next==face[f].nodeCount) next=0;
			patchArea=0.5*(faceNode(f,n)-centroid).cross(faceNode(f,next)-centroid);
			patchCentroid=1./3.*(faceNode(f,n)+faceNode(f,next)+centroid);
			face[f].centroid+=patchCentroid*patchArea.dot(face[f].normal);
			areaVec+=patchArea;
		}
		face[f].area=fabs(areaVec);
		face[f].centroid/=face[f].area;
	}
	
	cout << "[I Rank=" << Rank << "] Calculated face areas and centroids" << endl;
	
	// Loop through the cells and calculate the volumes
	double totalVolume=0.;
	for (int c=0;c<cellCount;++c) {
		double volume=0.;
		double patchVolume=0.;
		Vec3D height,basePatchArea,patchCentroid;
		int f,next;
		// Calculate cell centroid
		// First calculate an approximate one
		Vec3D centroid=0.;
		for (int cn=0;cn<cell[c].nodeCount;++cn) {
			centroid+=cellNode(c,cn);
		}
		centroid/=double(cell[c].nodeCount);
		// Break-up the volume into tetrahedras and add the volumes
		// Calculate the centroid of the cell by taking moments of each tetra
		cell[c].centroid=0.;
		double sign;
		for (int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			// Every cell face is broken to triangles
			for (int n=0;n<face[f].nodeCount;++n) {
				next=n+1;
				if (next==face[f].nodeCount) next=0;
				// Triangle area
				basePatchArea=0.5*(faceNode(f,n)-face[f].centroid).cross(faceNode(f,next)-face[f].centroid);
				// Height of the tetrahedra
				height=(face[f].centroid-centroid).dot(face[f].normal)*face[f].normal;
				// Fix face orientation issue
				sign=-1.;
				if (face[f].parent==c) sign=1.;
				patchVolume=sign*basePatchArea.dot(height)/3.;
				patchCentroid=0.25*(face[f].centroid+faceNode(f,n)+faceNode(f,next)+centroid);
				cell[c].centroid+=patchVolume*patchCentroid;
				// TODO Keep an eye on this
				//if (patchVolume<0.) cout << "[E Rank=" << Rank << "] Encountered error when calculating volume of cell " << c << endl;
				volume+=patchVolume;
			}
		}
		cell[c].volume=volume;
		cell[c].centroid/=volume;
		totalVolume+=volume;
	}
	
	cout << "[I Rank=" << Rank << "] Total Volume= " << setw(16) << setprecision(8) << scientific << totalVolume << endl;
	globalTotalVolume=0.;
	MPI_Allreduce (&totalVolume,&globalTotalVolume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if (Rank==0) cout << "[I] Global Total Volume= " << globalTotalVolume << endl;
	
	for (int f=0;f<faceCount;++f) {
		if (face[f].normal.dot(face[f].centroid-cell[face[f].parent].centroid)<=0.) {
			cout << "[W Rank=" << Rank << "] Face " << f << " normal is pointing in to its parent ... fixing " << endl;
			// Need to swap the face and reflect the area vector TODO Check if following is right
			face[f].normal*=-1.;
			vector<int>::reverse_iterator rit;
			face[f].nodes.assign(face[f].nodes.rbegin(),face[f].nodes.rend());
		}
	}
	
	return 0;

}


Node::Node(double x, double y, double z) {
	comp[0]=x;
	comp[1]=y;
	comp[2]=z;
}

Cell::Cell(void) {
	;
}

Node& Grid::cellNode(int c, int n) {
	return node[cell[c].nodes[n]];
};

Face& Grid::cellFace(int c, int f) {
	return face[cell[c].faces[f]];
};

Node& Grid::faceNode(int f, int n) {
	return node[face[f].nodes[n]];
};

void Grid::mpi_handshake(void) {
	
	sendCells.resize(np);
	recvCells.resize(np);
	vector<int> sendCount (np,0);
	vector<int> recvCount (np,0); 
	recvCount.resize(np);
	
	
	for (int g=0;g<ghostCount;++g) {
		int p=ghost[g].partition;
		recvCells[p].push_back(ghost[g].globalId);
	}
	for (int p=0;p<np;++p) recvCount[p]=recvCells[p].size();
	
	// I know which cells to receive from each processor
	// I don't know which cells to send
	// I don't even know how many to send
	// First communicate to figure out which processor request how many ghost cell data
	
	MPI_Alltoall(&recvCount[0],1,MPI_INT,&sendCount[0],1,MPI_INT,MPI_COMM_WORLD);
	

	for (int p=0;p<np;++p) {
		sendCells[p].resize(sendCount[p]);
		MPI_Sendrecv(&recvCells[p][0],recvCount[p],MPI_INT,p,0,
					 &sendCells[p][0],sendCount[p],MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}

	// Commit MPI_VEC3D
	MPI_Datatype types[2]={MPI_INT,MPI_DOUBLE};
	int block_lengths[2];
	MPI_Aint displacements[2];
	
	mpiGeomPack dummy4;
	displacements[0]=(long) &dummy4.ids[0] - (long) &dummy4;
	displacements[1]=(long) &dummy4.data[0] - (long) &dummy4;
	block_lengths[0]=2;
	block_lengths[1]=4;
	MPI_Type_create_struct(2,block_lengths,displacements,types,&MPI_GEOM_PACK);
	MPI_Type_commit(&MPI_GEOM_PACK);
	
}

void Grid::mpi_get_ghost_geometry(void) {
	
	for (int p=0;p<np;++p) {
		if (Rank!=p) {
			mpiGeomPack sendBuffer[sendCells[p].size()];
			mpiGeomPack recvBuffer[recvCells[p].size()];
			int id;
			for (int g=0;g<sendCells[p].size();++g) {
				id=maps.cellGlobal2Local[sendCells[p][g]];
				sendBuffer[g].ids[0]=myOffset+id;
				sendBuffer[g].ids[1]=id;
				for (int i=0;i<3;++i) sendBuffer[g].data[i]=cell[id].centroid[i];
				sendBuffer[g].data[3]=cell[id].volume;
			}
			
			MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GEOM_PACK,p,0,recvBuffer,recvCells[p].size(),MPI_GEOM_PACK,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			for (int g=0;g<recvCells[p].size();++g) {
				id=maps.ghostGlobal2Local[recvCells[p][g]];
				ghost[id].matrix_id=recvBuffer[g].ids[0];
				ghost[id].id_in_owner=recvBuffer[g].ids[1];
				for (int i=0;i<3;++i) ghost[id].centroid[i]=recvBuffer[g].data[i];
				ghost[g].volume=recvBuffer[g].data[3];
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	return;
} 
