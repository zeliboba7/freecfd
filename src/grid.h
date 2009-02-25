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
#include "vec3d.h"
using namespace std;
#include <parmetis.h>

class Node : public Vec3D {
public:
	unsigned int id, globalId;
	std::vector<int> cells;
	std::vector<int> ghosts;
	std::map<int,double> average;
	std::set<int> bcs;
	Node(double x=0., double y=0., double z=0.);
};

class Face {
public:
	int bc,parentIndex;
	unsigned int id;
	int parent,neighbor;
	unsigned int nodeCount;
	Vec3D centroid;
	Vec3D normal;
	std::map<int,double> average;
	double area;
	double mdot,weightL,weightR;
	std::vector<int> nodes;
	Node& node(int n);
};

class Cell {
public:
	unsigned int nodeCount,faceCount,neighborCellCount,ghostCount,globalId,globalCellCount;
	double volume,lengthScale,closest_wall_distance;
	Vec3D centroid;
	std::vector<int> nodes;
	std::vector<int> faces;
	std::vector<int> neighborCells;
	std::vector<int> ghosts;
	double p,T,rho;
	Vec3D v,grad[5];
	// Gradients are stored as p,u,v,w,T in order
	std::map<int,Vec3D> gradMap;
	double update[5];
	Cell(void);
	bool HaveNodes(unsigned int &nodelistsize, unsigned int nodelist[]) ;
	Node& node(int n);
	Face& face(int f);
};

class Ghost {
public:
	unsigned int partition;
	unsigned int globalId;
	unsigned int matrix_id;
	std::vector<unsigned int> cells;
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
	unsigned int nodeCount,cellCount,faceCount;
	unsigned int globalNodeCount,globalCellCount,globalFaceCount,ghostCount;
	std::vector<Node> node;
	std::vector<Face> face;
	std::vector<Cell> cell;
	std::vector<Ghost> ghost;
	std::vector<int> noSlipFaces;
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
	void gradients_turb();
	void limit_gradients_turb(void);
	void lengthScales(void);
};


class GridRawData { // The data will be destroyed after processing
public:
	std::vector<double> x,y,z;
	std::vector<int> cellConnIndex,cellConnectivity;
	std::vector< vector<int> > bocoConnIndex,bocoConnectivity;
	std::map<string,int> bocoNameMap;
};

class IndexMaps { // The data will be destroyed after processing
public:
	std::vector<unsigned int> cellOwner; // takes cell global id and returns the owner rank
	std::map<unsigned int,unsigned int> nodeGlobal2Local;
	std::map<unsigned int,unsigned int> cellGlobal2Local;
	std::map<unsigned int,unsigned int> ghostGlobal2Local;
	// This stores the bc region numbers that 'some' nodes touch to
	// Those 'some' nodes happen to be the first nodes in each bc face connectivity list
	std::map<unsigned int, vector<int> > nodeBCregions;
	// For those 'some' nodes, this stores the conn list of the corresponsing bc faces
	std::map<unsigned int, vector<set<unsigned int> > > nodeBCregionFaceConn;

	idxtype* adjIndex;
	idxtype* adjacency;
};


#endif
