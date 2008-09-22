/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include<cmath>
using namespace std;

#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;
extern int np, Rank;

extern GridRawData raw;

int Grid::CreateNodesCells() {
	// Create the nodes and cells for each partition
	// Run through the list of cells and check if it belongs to current partition
	// Loop through the cell's nodes
	// Mark the visited nodes so that no duplicates are created (nodeFound array).
	bool nodeFound[globalNodeCount]={false};
	unsigned int nodeMap[globalNodeCount]={-1};

	nodeCount=0;
	
	for (unsigned int c=0;c<globalCellCount;++c) {
		if (raw.cellMap[c]==Rank) {
			int cellNodeCount;
			if (c<globalCellCount-1) {
				cellNodeCount=raw.cellConnIndex[c+1]-raw.cellConnIndex[c];
			} else {
				elemNodeCount=raw.cellConnectivity.size()-raw.cellConnIndex[globalCellCount-1];
			}
			
			unsigned int cellNodes[cellNodeCount];
			for (unsigned int n=0;n<cellNodeCount;++n) {
				if (!nodeFound[raw.cellConnectivity[raw.cellConnIndex[c]+n]]) {
					Node temp;
					temp.id=nodeCount;
					temp.globalId=raw.cellConnectivity[raw.cellConnIndex[c]+n];
					temp.comp[0]=x[temp.globalId];
					temp.comp[1]=y[temp.globalId];
					temp.comp[2]=z[temp.globalId];
					node.push_back(temp);
					nodeFound[temp.globalId]=true;
					nodeMap[temp.globalId]=temp.id;
					++nodeCount;
				}
				cellNodes[n]=nodeMap[raw.cellConnectivity[raw.cellConnIndex[c]+n]];
			}

			Cell temp;
			temp.Construct(elemTypes[c],cellNodes);
			temp.globalId=c;
			cell.push_back(temp);

		}
	}


}

