/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

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
#include <sys/stat.h>
#include <sys/types.h>
#include "vec3d.h"
using namespace std;
#include <parmetis.h>

class Node : public Vec3D {
public:
	int id, globalId;
	std::vector<int> cells;
	std::vector<int> ghosts;
	std::map<int,double> average;
	std::set<int> bcs;
	Node(double x=0., double y=0., double z=0.);
};

class Face {
public:
	int bc,parentIndex;
	int id;
	int parent,neighbor;
	int nodeCount;
	Vec3D centroid;
	Vec3D normal;
	std::map<int,double> average;
	double area;
	double mdot,weightL;
	std::vector<int> nodes;
	Node& node(int n);
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
	double p,T,rho,dt;
	Vec3D v,grad[5];
	// Gradients are stored as p,u,v,w,T in order
	std::map<int,Vec3D> gradMap;
	double update[5];
	Cell(void);
	bool HaveNodes(int &nodelistsize, int nodelist[]) ;
	Node& node(int n);
	Face& face(int f);
};

class Ghost {
public:
	int partition;
	int globalId;
	int matrix_id;
	int id_in_owner;
	std::vector<int> cells;
	double p,T,rho,closest_wall_distance;
	// Gradients are stored as p,u,v,w,T in order
	Vec3D v,centroid,grad[5];
	double update[5]; // TODO do we need these updates for ghosts?
};

class Grid {
public:
	string fileName;
	int myOffset;
	vector<int> partitionOffset;
	int nodeCountOffset;
	int nodeCount,cellCount,faceCount;
	int globalNodeCount,globalCellCount,globalFaceCount,ghostCount;
	double globalTotalVolume;
	std::vector<int> boundaryFaceCount;
	std::vector<int> boundaryNodeCount;
	std::vector<int> globalBoundaryFaceCount;
	std::vector<int> globalBoundaryNodeCount;
	std::vector<Node> node;
	std::vector<Face> face;
	std::vector<Cell> cell;
	std::vector<Ghost> ghost;
	Grid();
	int read(string);
	int readCGNS();
	int readTEC();
	//int reorderRCM();
	int scale();
	int rotate();
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
	void gradMaps();
	void gradients();
 	void limit_gradients(void);
	void lengthScales(void);
};


class GridRawData { // The data will be destroyed after processing
public:
	std::vector<double> x,y,z;
	std::vector<int> cellConnIndex,cellConnectivity;
	std::vector< set<int> > bocoNodes; // Node list for each boundary condition region
	std::map<string,int> bocoNameMap;
};

class IndexMaps { // The data will be destroyed after processing
public:
	std::vector<int> cellOwner; // takes cell global id and returns the owner rank
	std::map<int,int> nodeGlobal2Local;
	std::map<int,int> cellGlobal2Local;
	std::map<int,int> ghostGlobal2Local;
	std::vector<int> nodeGlobal2Output;
	idxtype* adjIndex;
	idxtype* adjacency;
};


#endif
