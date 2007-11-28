#ifndef VECSPARSE_H
#define VECSPARSE_H

#include <vector>

class VecSparse {
public:
	std::vector<double> data;
	std::vector<int> indices;

	VecSparse(void);
	int position (int index);
	void put (int index, double a);
	double get(int index);
	void flush(void);
	VecSparse &operator*=(const double &);
	VecSparse operator*(const double &);
	VecSparse &operator/=(const double &);
	VecSparse operator/(const double &);
	VecSparse &operator+=(const double &);
	VecSparse operator+(const double &);
	VecSparse &operator-=(const double &);
	VecSparse operator-(const double &);

	VecSparse &operator+=(const VecSparse &);
	VecSparse &operator-=(const VecSparse &);
};

VecSparse operator*(const double &left, const VecSparse &right);
VecSparse operator/(const double &left, const VecSparse &right);
VecSparse operator+(const double &left, const VecSparse &right);
VecSparse operator-(const double &left, const VecSparse &right);

#endif
