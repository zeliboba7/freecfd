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
#include "polynomial.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

void Polynomial::set(string type_in,int piece_count,vector<double> list) {
	// Takes in a list of numbers containing beginning range and coefficients for each piece
	begin.resize(piece_count);
	coeff.resize(piece_count);
	int order=list.size()/piece_count-2; // e.g. a 2nd order poly has 4 entries (begin,coeff1,coeff2,coeff3)
	if (type_in=="schomate") type=SCHOMATE;
	else type=REGULAR;
	for (int p=0;p<piece_count;++p) {
		coeff[p].resize(order+1);
		int index=p*(order+2);
		begin[p]=list[index];
		for (int c=0;c<coeff[p].size();++c) coeff[p][c]=list[index+1+c];
	}
	
	return;
}

double Polynomial::eval(double x) {
	double sum=0.;
	int p;
	
	if (type==REGULAR) {
		for (p=0;p<begin.size()-1;p++) {
			if (x<begin[p+1]) {
				// for each coefficient
				sum+=coeff[p][0];
				for (int c=1;c<coeff.size();++c) {
					sum+=coeff[p][c]*pow(x,double(c));	
				}
				return sum;
			}
		}
		// If function didn't return above, x is at the last piece
		// for each coefficient
		p=begin.size()-1;
		sum+=coeff[p][0];
		for (int c=1;c<coeff.size();++c) {
			sum+=coeff[p][c]*pow(x,double(c));	
		}
		return sum;
	} else if (type==SCHOMATE) {
		for (p=0;p<begin.size()-1;p++) {
			if (x<begin[p+1]) {
				x/=1000.;
				return coeff[p][0]+coeff[p][1]*x+coeff[p][2]*x*x+coeff[p][3]*x*x*x+coeff[p][4]/(x*x);
			}
		}
		// If function didn't return above, x is at the last piece
		p=begin.size()-1;
		x/=1000.;
		return coeff[p][0]+coeff[p][1]*x+coeff[p][2]*x*x+coeff[p][3]*x*x*x+coeff[p][4]/(x*x);
	}
}
