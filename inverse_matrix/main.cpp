#include <iostream>
#include <iomanip>
#include <ctime>
#include "funcs.h"
using namespace std;

int main()
{
	clock_t start_time, end_time;
	double time;
	int ndim;
	cout << "输入矩阵的维度(10以上不显示矩阵):" << flush;
	cin >> ndim;
	Matrix A(0,1,ndim);
	A.creat();
	cout << "原矩阵:" << endl;
	if (ndim <= 10) A.display(0); 
	cout << endl;

	start_time = clock();
	A.GJ_Inverse();
	end_time = clock();
	time = double(end_time - start_time);
	cout << "高斯消元法,用时" << setw(5) << time/CLOCKS_PER_SEC << "s,得出逆矩阵:" << endl;
	if (ndim <= 10) A.display(1);
	cout << endl;

	start_time = clock();
	A.LU_Inverse();
	end_time = clock();
	time = double(end_time - start_time);
	cout << "LU分解法,用时" << setw(5) << time/CLOCKS_PER_SEC << "s,得出逆矩阵:" << endl;
	if (ndim <= 10) A.display(1);
	cout << endl;

	system("pause");
	return 0;
}