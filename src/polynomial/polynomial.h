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
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cmath>
#include <string>
using namespace std;

// Options for polynomial types
#define REGULAR 1
#define SCHOMATE 2

class Polynomial {
public:
	// A piece-wise polynomial
	int type;
	// Beginning point of each piece
	vector<double> begin;
	// Coefficients for each piece
	vector<vector<double> > coeff;
	void set(string type_in,int piece_count,vector<double> list);
	double eval(double x);
};

#endif
