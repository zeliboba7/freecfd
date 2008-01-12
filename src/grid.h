#ifndef GRID_H
#define GRID_H

/*
  Classes for reading and storing grid information
*/

#include <cgnslib.h>
#include <vector>
#include <mpi.h>
#include "vec3d.h"
#include "sparse.h"


class Node : public Vec3D {
public:
	unsigned int id, globalId;
	std::vector<int> cells;
	Node(double x=0., double y=0., double z=0.);
};

class Face {
public:
	int bc;
	unsigned int id;
	int parent,neighbor;
	unsigned int nodeCount;
	Vec3D centroid;
	Vec3D normal;
	VecSparse cellContributions;
	double area;
	std::vector<int> nodes;
	Node& node(int n);
};

class Cell {
public:
	unsigned int nodeCount,faceCount,globalId,globalCellCount;
	ElementType_t type;
	double volume;
	Vec3D centroid;
	std::vector<int> nodes;
	std::vector<int> faces;
	double rho;
	Vec3D v,grad[5];
	double p;
	double flux[5];
	Cell(void);
	int Construct(const ElementType_t elemType,unsigned int nodeList[]);
	bool HaveNodes(unsigned int &nodelistsize, unsigned int nodelist[]) ;
	Node& node(int n);
	Face& face(int f);
};

class Ghost {
public:
	unsigned int partition;
	unsigned int globalId;
	double rho;
	Vec3D v;
	double p;
};

class Grid {
public:
	string fileName;
	unsigned int nodeCount,cellCount,faceCount;
	unsigned int globalNodeCount,globalCellCount,globalFaceCount,ghostCount;
	std::vector<Node> node;
	std::vector<Face> face;
	std::vector<Cell> cell;
	std::vector<Ghost> ghost;
	int ReadCGNS();
	Grid();
	int read(string);
	int face_exists(int &parentCell);
	void gradients();
};

#endif
