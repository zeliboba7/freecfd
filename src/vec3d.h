#ifndef VEC3D_H
#define VEC3D_H

#include <iostream>
using namespace std;

class Vec3D {
public:
  double comp[3];
	Vec3D(double x=0., double y=0., double z=0.);
	double dot(const Vec3D &right);
	Vec3D cross(const Vec3D &right);
	Vec3D &operator=(const Vec3D &);
	Vec3D &operator=(const double &);
	Vec3D &operator*=(const double &);
	Vec3D operator*(const double &);
	Vec3D &operator/=(const double &);
	Vec3D operator/(const double &);
	Vec3D &operator+=(const double &);
	Vec3D &operator+=(const Vec3D &);
	Vec3D operator+(const double &);
	Vec3D &operator-=(const double &);
	Vec3D &operator-=(const Vec3D &);
	Vec3D operator-(const double &);
	bool operator==(const Vec3D &);
	bool operator!=(const Vec3D &);
};

double fabs(const Vec3D vec);
Vec3D operator*(const double &left, const Vec3D &right);
Vec3D operator/(const double &left, const Vec3D &right);
Vec3D operator+(const double &left, const Vec3D &right);
Vec3D operator-(const double &left, const Vec3D &right);

Vec3D operator+(const Vec3D &left, const Vec3D &right);
Vec3D operator-(const Vec3D &left, const Vec3D &right);

ostream &operator<<(ostream &output,const Vec3D &right);

#endif
