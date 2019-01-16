#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
using namespace std;

int accept = 0;
int total = 0;

int main()
{
	double fRand(double, double);
	double *next_point(double *, int, double, double, string method = "period");
	double func_integrate(int, double);

	//Markov链系数
	int D = 5;
	double a = 2.0;
	double *x;
	int i, j;
	extern int accept;
	extern int total;
	

	////游走,然后输出到屏幕和文件
	//ofstream outfile("data1.csv", ios::out);
	////cout << setiosflags(ios::fixed) << setiosflags(ios::left) << setprecision(4);
	//for (j = 0; j < d; j++) {
	//	//cout << setw(8) << x[j];
	//	if (j < d - 1) outfile << x[j] << ",";
	//	else outfile << x[j] << endl;
	//}
	////cout << endl;
	//for (i = 0; i < 20000; i++) {
	//	x = next_point(x, d, a, 0.1*a, "period");
	//	for (j = 0; j < d; j++) {
	//		//cout << setw(8) << x[j];
	//		if (j < d - 1) outfile << x[j] << ",";
	//		else outfile << x[j] << endl;
	//	}
	//	//cout << endl;
	//}
	//outfile.close();

	////(2)N_equi=10000, N0=50, m=10000
	////精确值
	//cout << "精确值:" << func_integrate(5, 2) << endl;

	//int L = 400; //运行次数
	//ofstream outfile("data2.csv", ios::app);
	//for (int k = L; k < L + 15; k++) {
	//	srand(k);
	//	
	//	//新建随机值x[D]
	//	x = new double[D];
	//	for (i = 0; i < D; i++)
	//		*(x + i) = fRand(0, a);

	//	int iter_sum = 0, iter = 0;
	//	double sum = 0, prod, result;

	//	
	//	while (iter_sum < 10000) {
	//		x = next_point(x, 5, 2, 0.2, "period");
	//		iter += 1;
	//		if (iter < 10000) continue;
	//		if ((iter - 10000) % 50 == 0) {
	//			prod = 1;
	//			for (j = 0; j < 5; j++) prod *= (1 + x[j] / 2);

	//			sum += prod;
	//			iter_sum += 1;
	//		}
	//	}
	//	result = pow((1 - exp(-2)), 5)*sum / 10000;
	//	cout << "计算值:" << result << endl;
	//	outfile << k << ',' << result << endl;
	//	delete[] x;
	//	x = NULL;
	//}
	//outfile.close();

	//(4)N_equi=2000~16000, N0=50, m=1e4
	//精确值
	cout << "精确值:" << func_integrate(D, 2) << endl;

	int L = 160; //运行次数

	ofstream outfile("data6_0.05a.csv", ios::app);
	for (int k = L; k < L + 40; k++) {
		srand(k + 4600);

		//新建随机值x[D]
		x = new double[D];
		for (i = 0; i < D; i++)
			*(x + i) = fRand(0, a);

		int iter_sum = 0, iter = 0;
		double sum = 0, prod, result;


		while (iter_sum < 10000) {
			x = next_point(x, D, 2, 0.1, "period");
			iter += 1;
			if (iter < 10000) continue;
			if ((iter - 10000) % 50 == 0) {
				prod = 1;
				for (j = 0; j < D; j++) prod *= (1 + x[j] / 2);

				sum += prod;
				iter_sum += 1;
			}
		}
		result = pow((1 - exp(-2)), D)*sum / 10000;
		cout << "计算值(" << k << "): " << result << ',' << double(accept) / double(total) << endl;
		outfile << k << ',' << result << ',' << double(accept) / double(total) << endl;
		delete[] x;
		x = NULL;
	}
	outfile.close();

}

double fRand(double fMin, double fMax)
{
	double t = double(rand()) / (double(RAND_MAX) + 1.0);
	return fMin + t * (fMax - fMin);
}

double *next_point(double *x0, int D, double a, double delta, string method)
{
	//给出马尔可夫链的下一个点
	//x0:输入; D:维度; a:取值上界; delta:delta; method:"reflect" or "period";
	double fRand(double, double);
	int i;
	double *x_temp, *dx, e_sum;
	extern int accept;
	extern int total;

	x_temp = new double[D];
	dx = new double[D];

	for (i = 0; i < D; i++) {
		dx[i] = delta * (fRand(0, 1) - 0.5) * 2;
		x_temp[i] = x0[i] + dx[i];
	}

	//反射折回or周期折回
	if (method == "reflect") {
		for (i = 0; i < D; i++) {
			if (x_temp[i] > a) x_temp[i] = 2 * a - x_temp[i];
			else if (x_temp[i] < 0) x_temp[i] = -x_temp[i];
		}
	}
	else {
		for (i = 0; i < D; i++) {
			if (x_temp[i] > a) x_temp[i] -= a;
			else if (x_temp[i] < 0) x_temp[i] += a;
		}
	}
	
	//是否更新值
	total += 1;
	e_sum = 0;
	for (i = 0; i < D; i++) {
		e_sum += x_temp[i] - x0[i];
	}
	if (exp(-e_sum) >= fRand(0, 1)){
		accept += 1;
		for (i = 0; i < D; i++)
			x0[i] = x_temp[i];
	}

	delete x_temp, dx;
	x_temp = NULL;
	dx = NULL;
	return x0;
}

double func_integrate(int D, double a)
{
	return pow((1.5 - (1.5 + a / 2)*exp(-a)), D);
}