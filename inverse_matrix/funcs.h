#include <stdlib.h>

class Matrix
{
public:
	Matrix() { fMin = 0.0; fMax = 1.0; n = 5; }
	Matrix(double min, double max, int dim) :fMin(min), fMax(max), n(dim) {}
	void creat();
	void display(int op = 0);
	void GJ_Inverse();
	void LU_Inverse();
private:
	double fRand();
	double **M, **M_raw, **eye;
	double fMin, fMax;
	int n;
};