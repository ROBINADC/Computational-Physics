/*
This head includes fRand, creat_m for creating matrix, creat_v for creating vector,
GJElimination, print_m for printing matrix, print_v for printing vector, copy_m for copying matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fRand(double fMin, double fMax)
{
	double t = (double)rand() / ((double)RAND_MAX + 1.0);
	return fMin + t * (fMax - fMin);
}

double **creat_m(int n, double fMin, double fMax)
{
	//creat a matrix that includes coefficient terms and constant terms, like N * (N + 1)
	//Also it is filled up with random numbers range from fMin to fMax.
	double **A;
	int i, j;
	A = (double **)malloc(sizeof(double *)*n);
	for (i = 0; i < n; i++) {
		*(A+i) = (double *)malloc(sizeof(double)*(n + 1));
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j <= n; j++) {
			A[i][j] = fRand(fMin, fMax);
		}
	}
	return A;
}

double *creat_v(int n, double fMin, double fMax)
{
	//it seems that this is useless
	return 0;
}

void print_m(double **A, int n)
{
	//print a matrix with N * (N + 1)
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j <= n; j++) {
			if (j == n)
				printf("| ");
			printf("%10.6f ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_v(double *vec, int n)
{
	//print a vector with dimension n
	for (int i = 0; i < n; i++) {
		printf("%10.6f\n", vec[i]);
	}
	printf("\n");
}

double **copy_m(double **A, int n)
{
	//copy an N * (N + 1) matrix
	double **A_copy;
	A_copy = (double **)malloc(sizeof(double *)*n);
	for (int i = 0; i < n; i++)
		*(A_copy + i) = (double *)malloc(sizeof(double)*(n + 1));
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= n; j++)
			A_copy[i][j] = A[i][j];
	return A_copy;
}

struct result
{
	double **AT;	//三角矩阵
	double *x;	//解
	double *d;	//解的偏离
};

result GJElimination(double **A, int n)
{
	//make tri-matrix, find solution x, calculate bias c, and return result.
	result rst;
	int i, j, k, n_max;
	double max, swap, m, temp;

	//初始化结构体
	rst.AT = (double **)malloc(sizeof(double *)*n);
	for (i = 0; i < n; i++) {
		*(rst.AT + i) = (double *)malloc(sizeof(double)*(n + 1));
	}
	rst.x = (double *)malloc(sizeof(double)*n);
	rst.d = (double *)malloc(sizeof(double)*n);

	rst.AT = copy_m(A, n);
	for (j = 0; j < n - 1; j++) {	//j次数，执行行数
		max = fabs(rst.AT[j][j]);
		n_max = j;
		for (i = j; i < n; i++) {	//按当前列找到最大行
			if (fabs(rst.AT[i][j]) > max) {
				max = fabs(rst.AT[i][j]);
				n_max = i;
			}
		}

		if (fabs(max) < 1.0e-6) {	//是否为奇异矩阵
			printf("sigular matrix: %lf\n", fabs(max));
		}

		for (k = j; k <= n; k++) {	//k一般是列标
			//交换最大行与当前行，由于列标比j小的均为0，仅交换>=j的即可
			swap = rst.AT[j][k];
			rst.AT[j][k] = rst.AT[n_max][k];
			rst.AT[n_max][k] = swap;
		}

		for (i = j + 1; i < n; i++) {	//j的下一行直至末行
			m = rst.AT[i][j] / rst.AT[j][j];	//构建乘子
			for (k = j; k <= n; k++) {
				rst.AT[i][k] -= m * rst.AT[j][k];
			}
		}
	}

	//seek for result
	rst.x[n - 1] = rst.AT[n - 1][n] / rst.AT[n - 1][n - 1];
	for (i = n - 2; i >= 0; i--) {
		temp = rst.AT[i][n];
		for (k = i + 1; k < n; k++) {
			temp -= rst.AT[i][k] * rst.x[k];
		}
		rst.x[i] = temp / rst.AT[i][i];
	}

	//计算解的误差
	for (i = 0; i < n; i++) {
		rst.d[i] = rst.AT[i][n];
		for (j = 0; j < n; j++) {
			rst.d[i] -= rst.AT[i][j] * rst.x[j];	//用结构体时刻注意用成员
		}
	}

	return rst;
}