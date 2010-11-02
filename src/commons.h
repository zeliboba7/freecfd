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
#ifndef COMMONS_H
#define COMMONS_H

#define NONE -1
// Equation options
#define NS 1
#define HEAT 2
// Gradient test options
#define LINEAR 1
#define QUADRATIC 2

extern int Rank,np;
extern int gradient_test;
extern double min_x,max_x;

#endif
