#ifndef BC_H
#define BC_H

#include <vector>

class BCregion {
public:
	string type;
	double rho,p;
	Vec3D v;
};

class BC {
public:
	vector<BCregion> region;
};

#endif
