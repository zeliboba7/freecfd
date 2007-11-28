#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
#include "vec3d.h"

Vec3D::Vec3D (double x, double y, double z) {
  comp[0]=x;
  comp[1]=y;
  comp[2]=z;
}

double Vec3D::dot(const Vec3D &right) {
 return (comp[0]*right.comp[0]+comp[1]*right.comp[1]+comp[2]*right.comp[2]);
}

double fabs(const Vec3D vec) {
 return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]+vec.comp[2]*vec.comp[2]);
}

Vec3D Vec3D::cross(const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=comp[1]*right.comp[2]-comp[2]*right.comp[1];
 temp.comp[1]=-comp[0]*right.comp[2]+comp[2]*right.comp[0];
 temp.comp[2]=comp[0]*right.comp[1]-comp[1]*right.comp[0]; 
 return temp;
}

Vec3D &Vec3D::operator=(const Vec3D &right) {
 comp[0]=right.comp[0];
 comp[1]=right.comp[1];
 comp[2]=right.comp[2];
return *this;
}

Vec3D &Vec3D::operator=(const double &right) {
 comp[0]=right;
 comp[1]=right;
 comp[2]=right;
return *this;
}

Vec3D &Vec3D::operator*=(const double &right) {
 comp[0]*=right;
 comp[1]*=right;
 comp[2]*=right;
return *this;
}

Vec3D Vec3D::operator*(const double &right) {
 Vec3D temp;
 temp.comp[0]=comp[0]*right;
 temp.comp[1]=comp[1]*right;
 temp.comp[2]=comp[2]*right;
return temp;
}

Vec3D &Vec3D::operator/=(const double &right) {
 comp[0]/=right;
 comp[1]/=right;
 comp[2]/=right;
return *this;
}

Vec3D Vec3D::operator/(const double &right) {
 Vec3D temp;
 temp.comp[0]=comp[0]/right;
 temp.comp[1]=comp[1]/right;
 temp.comp[2]=comp[2]/right;
return temp;
}

Vec3D &Vec3D::operator+=(const double &right) {
 comp[0]+=right;
 comp[1]+=right;
 comp[2]+=right;
return *this;
}

Vec3D &Vec3D::operator+=(const Vec3D &right) {
 comp[0]+=right.comp[0];
 comp[1]+=right.comp[1];
 comp[2]+=right.comp[2];
return *this;
}

Vec3D Vec3D::operator+(const double &right) {
 Vec3D temp;
 temp.comp[0]=comp[0]+right;
 temp.comp[1]=comp[1]+right;
 temp.comp[2]=comp[2]+right;
return temp;
}

Vec3D &Vec3D::operator-=(const double &right) {
 comp[0]-=right;
 comp[1]-=right;
 comp[2]-=right;
return *this;
}

Vec3D &Vec3D::operator-=(const Vec3D &right) {
 comp[0]-=right.comp[0];
 comp[1]-=right.comp[1];
 comp[2]-=right.comp[2];
return *this;
}

Vec3D Vec3D::operator-(const double &right) {
 Vec3D temp;
 temp.comp[0]=comp[0]-right;
 temp.comp[1]=comp[1]-right;
 temp.comp[2]=comp[2]-right;
return temp;
}

Vec3D operator*(const double &left, const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=left*right.comp[0];
 temp.comp[1]=left*right.comp[1];
 temp.comp[2]=left*right.comp[2];
return temp;
}

Vec3D operator/(const double &left, const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=left/right.comp[0];
 temp.comp[1]=left/right.comp[1];
 temp.comp[2]=left/right.comp[2];
return temp;
}

Vec3D operator+(const double &left, const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=left+right.comp[0];
 temp.comp[1]=left+right.comp[1];
 temp.comp[2]=left+right.comp[2];
return temp;
}

Vec3D operator+(const Vec3D &left, const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=left.comp[0]+right.comp[0];
 temp.comp[1]=left.comp[1]+right.comp[1];
 temp.comp[2]=left.comp[2]+right.comp[2];
return temp;
}

Vec3D operator-(const double &left, const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=left-right.comp[0];
 temp.comp[1]=left-right.comp[1];
 temp.comp[2]=left-right.comp[2];
return temp;
}

Vec3D operator-(const Vec3D &left, const Vec3D &right) {
 Vec3D temp;
 temp.comp[0]=left.comp[0]-right.comp[0];
 temp.comp[1]=left.comp[1]-right.comp[1];
 temp.comp[2]=left.comp[2]-right.comp[2];
return temp;
}

bool Vec3D::operator==(const Vec3D &right) {
	return (comp[0]==right.comp[0] && comp[1]==right.comp[1] && comp[2]==right.comp[2] );
}

bool Vec3D::operator!=(const Vec3D &right) {
	return (comp[0]!=right.comp[0] | comp[1]!=right.comp[1] | comp[2]!=right.comp[2] );
}

ostream &operator<<(ostream &output,const Vec3D &right) {
  output << "{" << right.comp[0] << "," << right.comp[1] << "," << right.comp[2] << "}";
  return output;
}
