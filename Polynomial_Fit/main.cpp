#include <iostream>
#include <cmath>
#include <iomanip>
#include "funcs.h"
using namespace std;

int main()
{
	double x[11] = { 1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000 };
	double y[11] = { 75.995, 91.972,105.711, 123.203,131.669,150.697,179.323,203.212,226.505,
				249.633,281.422 };
	double *s, *t, **L;
	int i, j, k;
	result p;

	s = new double[2 * 3 + 1];
	t = new double[3 + 1];
	L = new double *[3 + 1];
	for (i = 0; i < 3 + 1; i++) {
		L[i] = new double[5];
	}
	
	//赋值给s
	for (k = 0; k < 7; k++) {
		s[k] = 0;
		for (i = 0; i < 11; i++) {
			s[k] += pow(x[i], k);
		}
	}

	//赋值给t
	for (k = 0; k < 4; k++) {
		t[k] = 0;
		for (i = 0; i < 11; i++) {
			t[k] += y[i] * pow(x[i], k);
		}
	}

	//形成增广矩阵L
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 5; j++) {
			if (j < 4) L[i][j] = s[i + j];
			else if (j = 4) L[i][j] = t[i];
		}
	}

	p = GJElimination(L, 4);
	cout << "1900到2000年的人口数:" << endl;
	cout << "年份" << "     " << "人口数(百万)" << endl;
	for (i = 0; i < 11; i++) {
		cout << x[i] << "     " << y[i] << endl;
	}
	cout << "\n拟合曲线:(y:人口数(百万);x:年份)" << endl;
	cout << "y = " << p.x[0] << " + " << p.x[1] << "x + " << p.x[2] << "x^2 + " << p.x[3] << "x^3" << endl;

	double pred;
	pred = p.x[0] + p.x[1] * 2010 + p.x[2] * 2010 * 2010 + p.x[3] * 2010 * 2010 * 2010;
	cout << "\n以此预测美国2010年的人口数为:" << pred <<"百万"<< endl;
	cout << "美国2010年的实际人口数:309.35百万" << endl;
	system("pause");
}