#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "funcs.h"

int main(int argc, char *argv[]) {
	double **A;
	int n;
	double fMin, fMax;
	result R;

	printf("��������ά��:");
	scanf_s("%d", &n);
	printf("\n����ϵ����Сֵ�����ֵ(�ո�ָ�):");
	scanf_s("%lf%lf", &fMin, &fMax);
	
	A = creat_m(n, fMin, fMax);
	printf("\n�Ѵ�������:\n");
	print_m(A, n);

	R = GJElimination(A, n);
	printf("������任Ϊ�����Ǿ���:\n");
	print_m(R.AT, n);
	printf("��:\n");
	print_v(R.x, n);
	printf("���ƫ��:\n");
	print_v(R.d, n);

	system("pause");
	return 0;
}