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
	cout << "��������ά��(10���ϲ���ʾ����):" << flush;
	cin >> ndim;
	Matrix A(0,1,ndim);
	A.creat();
	cout << "ԭ����:" << endl;
	if (ndim <= 10) A.display(0); 
	cout << endl;

	start_time = clock();
	A.GJ_Inverse();
	end_time = clock();
	time = double(end_time - start_time);
	cout << "��˹��Ԫ��,��ʱ" << setw(5) << time/CLOCKS_PER_SEC << "s,�ó������:" << endl;
	if (ndim <= 10) A.display(1);
	cout << endl;

	start_time = clock();
	A.LU_Inverse();
	end_time = clock();
	time = double(end_time - start_time);
	cout << "LU�ֽⷨ,��ʱ" << setw(5) << time/CLOCKS_PER_SEC << "s,�ó������:" << endl;
	if (ndim <= 10) A.display(1);
	cout << endl;

	system("pause");
	return 0;
}