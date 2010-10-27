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
#ifndef GRID_H
#define GRID_H

// Face bc type numbering
#define INTERNAL_FACE -1
#define UNASSIGNED_FACE -2
#define GHOST_FACE -3
/*
  Classes for reading and storing grid information
*/

#include <cgnslib.h>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <mpi.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;
#include <parmetis.h>

#include "vec3d.h"

class GridRawData { // This data will be destroyed after processing
public:
	std::vector<Vec3D> node;
	std::vector<int> cellConnIndex,cellConnectivity;
	std::vector< set<int> > bocoNodes; // Node list for each boundary condition region
	std::map<string,int> bocoNameMap;
};

class IndexMaps { // This data will be destroyed after processing
public:
	std::vector<int> cellOwner; // takes cell global id and returns the owner rank
	std::map<int,int> nodeGlobal2Local;
	std::map<int,int> cellGlobal2Local;
	std::map<int,int> ghostGlobal2Local;
	std::vector<int> nodeGlobal2Output;
	std::vector<int> face2bc; // face index to bc array index map
	idxtype* adjIndex;
	idxtype* adjacency;
};

class Node : public Vec3D {
public:
	int id, globalId; // id is the local index in the current processor
	std::vector<int> cells; // list of cells (ids) sharing this node
	std::vector<int> ghosts; // list of ghosts (inter-partition cells) sharing this node 
	std::vector<int> faces; // list of faces touching this node
	std::map<int,double> average; // indices of cells in the averaging stencil and corresponding weights
	Node(double x=0., double y=0., double z=0.);
};

class Face {
public:
	int bc; // A face can only be on one bc zone
	int id;
	int parent,neighbor; // parent is the id of the cell owning this face
	// A cell owns a face if the face normal is pointing outwards from it
	// The other cell is called the neighbor
	int nodeCount;
	Vec3D centroid;
	Vec3D normal; // This should point outwards from the parent cell center
	std::map<int,double> average; // indices of cells in the averaging stenceil and corresponding weights
	double area; 
	std::vector<int> nodes; // Nodes of this face
};

class Cell {
public:
	int nodeCount,faceCount,neighborCellCount,ghostCount,globalId,globalCellCount;
	double volume,lengthScale,closest_wall_distance;
	Vec3D centroid;
	std::vector<int> nodes;
	std::vector<int> faces;
	std::vector<int> neighborCells;
	std::vector<int> ghosts;
	Cell(void);
	bool HaveNodes(int &nodelistsize, int nodelist[]) ;
};

class Ghost {
public:
	int partition;
	int globalId;
	int matrix_id;
	int id_in_owner;
	std::vector<int> cells;
	Vec3D centroid;
	double volume;
};

class Grid {
public:
	int gid;
	int dimension; // 2 or 3
	int bcCount;
	double lengthScale;
	GridRawData raw;
	IndexMaps maps;
	string fileName;
	int myOffset,Rank,np;
	vector<int> partitionOffset;
	int nodeCountOffset;
	int nodeCount,cellCount,faceCount;
	int globalNodeCount,globalCellCount,globalFaceCount,ghostCount;
	double globalTotalVolume;
	std::vector<std::vector<int> > boundaryFaceCount; // for each bc region in each proc 
	std::vector<int> globalBoundaryFaceCount;
	std::vector<Node> node;
	std::vector<Face> face;
	std::vector<Cell> cell;
	std::vector<Ghost> ghost;
	// Maps for MPI exchanges
	std::vector< std::vector<int> > sendCells;
	std::vector< std::vector<int> > recvCells;
	MPI_Datatype MPI_GEOM_PACK;
	Grid();
	void read(string);
	void setup(void);
	int readCGNS();
	int translate(Vec3D begin, Vec3D end);
	int scale(Vec3D anchor, Vec3D factor);
	int rotate(Vec3D anchor, Vec3D axis, double angle);
	int partition();
	int mesh2dual();
	int create_nodes_cells();
	int create_faces();
	int create_ghosts();
	void trim_memory();
	int areas_volumes();
	void nodeAverages();
	void sortStencil(Node& n);
	void interpolate_tetra(Node& n);
	void interpolate_tri(Node& n);
	void interpolate_line(Node& n);
	void nodeAverages_idw();
	void faceAverages();
	void mpi_handshake(void);
	void mpi_get_ghost_geometry(void);
	void write_cas(void);

	Node& cellNode(int c, int n);
	Face& cellFace(int c, int f);
	Node& faceNode(int f, int n);
};

// Custom MPI type to exhange ghost centroids
struct mpiGeomPack {
	int ids[2]; // contains globalId and matrix_id;
	double data[4];
};

#endif
