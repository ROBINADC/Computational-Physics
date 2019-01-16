#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "funcs.h"

int main(int argc, char *argv[]) {
	double **A;
	int n;
	double fMin, fMax;
	result R;

	printf("输入矩阵的维度:");
	scanf_s("%d", &n);
	printf("\n输入系数最小值和最大值(空格分隔):");
	scanf_s("%lf%lf", &fMin, &fMax);
	
	A = creat_m(n, fMin, fMax);
	printf("\n已创建矩阵:\n");
	print_m(A, n);

	R = GJElimination(A, n);
	printf("将矩阵变换为上三角矩阵:\n");
	print_m(R.AT, n);
	printf("解:\n");
	print_v(R.x, n);
	printf("解的偏差:\n");
	print_v(R.d, n);

	system("pause");
	return 0;
}