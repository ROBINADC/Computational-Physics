#include <iostream>
#include <iomanip>
#include <stdlib.h>
using namespace std;

class Matrix
{
public:
	Matrix() { fMin = 0.0; fMax = 1.0; n = 5; }
	Matrix(double min, double max, int dim) :fMin(min), fMax(max), n(dim) {}
	void creat();
	void display(int op=0);
	void GJ_Inverse();
	void LU_Inverse();
private:
	double fRand();
	double **M, **M_raw, **eye;	//MΪ������ľ���,M_rawΪԭ����,eyeΪ��λ����
	double fMin, fMax;
	int n;
};

double Matrix::fRand()
{
	return fMin + (double(rand()) / (double(RAND_MAX) + 1))*(fMax - fMin);
}

void Matrix::creat()
{
	//creat an n*n matrix M_raw, copy it to M, initialize eye
	int i, j;
	M = new double *[n];
	for (i = 0; i < n; i++) {
		M[i] = new double[n];
	}
	M_raw = new double *[n];
	for (i = 0; i < n; i++) {
		M_raw[i] = new double[n];
	}
	eye = new double *[n];
	for (i = 0; i < n; i++) {
		eye[i] = new double[n];
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			M_raw[i][j] = fRand();
			M[i][j] = M_raw[i][j];
		}
	}
}

void Matrix::display(int op)
{
	//print matrix: if op=0, print raw matrix; else, print modified matrix
	int i, j;
	cout << setiosflags(ios::fixed) << setiosflags(ios::right) << setprecision(4);

	if (op == 0) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << setw(8) << M_raw[i][j];
			}
			cout << endl;
		}
	}
	else if (op) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << setw(8) << M[i][j];
			}
			cout << endl;
		}
	}
}

void Matrix::GJ_Inverse()
{
	//complete inverse an matrix by Gauss-Jordan Elimination, and the result would be put in M
	int i, j, k, n_max;
	double max, swap, swap1, mul;

	//��M_raw��ʼ��M
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			M[i][j] = M_raw[i][j];

	//creat identity matrix
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) 
			eye[i][j] = i == j ? 1 : 0;
	}

	//�������Ϊ�����Ǿ���
	for (k = 0; k < n - 1; k++) {
		max = fabs(M[k][k]);
		n_max = k;
		for (i = k; i < n; i++) {
			if (fabs(M[i][k]) > max) {
				max = fabs(M[i][k]);
				n_max = i;
			}
		}

		for (j = 0; j < n; j++) {
			//������
			swap = M[k][j];
			swap1 = eye[k][j];
			M[k][j] = M[n_max][j];
			eye[k][j] = eye[n_max][j];
			M[n_max][j] = swap;
			eye[n_max][j] = swap1;
		}

		for (i = k + 1; i < n; i++) {
			//��Ϊ�����Ǿ���
			mul = M[i][k] / M[k][k];
			for (j = 0; j < n; j++) {
				M[i][j] -= M[k][j] * mul;
				eye[i][j] -= eye[k][j] * mul;
			}
		}
	}

	for (i = n - 1; i >= 0; i--) {	//�����һ�п�ʼ,ֱ����һ��
		for (j = 0; j < n; j++) {
			//�������Ǿ���ĶԽ�Ԫ��Ϊ1
			//M[i][j] /= M[i][i] ��ʱĬ���ڸ��еĶԽ�ԪΪ1��������Ҫִ��
			eye[i][j] /= M[i][i];
		}

		//��ԭ����Ϊ��λ����(ֻ������������Ҳ����)
		if (i > 0) {	//�����һ�п�ʼ�����ڶ���
			for (k = i; k < n; k++) {
				mul = M[i - 1][k];	//���ڶԽ�ԪΪ1���˼�����
				for (j = 0; j < n; j++) {
					eye[i - 1][j] -= mul * eye[k][j];
				}
			}
		}
	}

	//��������Ƶ�M
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			M[i][j] = eye[i][j];
}

void Matrix::LU_Inverse()
{
	//LU�ֽ���������,����һ������������ԭ����M
	double **L, **U, **L_, **U_;
	int i, j, k;
	double t = 1.0; //�ж�������Ƿ����

	//��M_raw��ʼ��M
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			M[i][j] = M_raw[i][j];

	/*
	for (j = 0; j < n; j++) {  //�������ĵ�һ��
		M[0][j] = M[0][j];
	}
	*/
	for (i = 1; i < n; i++) {
		M[i][0] = M[i][0] / M[0][0];	//�������ĵ�һ��
	}
	
	for (k = 1; k < n; k++) {
		for (j = k; j < n; j++) {
			double sum = 0.0;
			for (i = 0; i < k; i++) {
				sum += M[k][i] * M[i][j];
			}
			M[k][j] -= sum; //����U����Ĳ���
		}
		for (i = k + 1; i < n; i++) {
			double sum = 0.0;
			for (j = 0; j < k; j++) {
				sum += M[i][j] * M[j][k];
			}
			M[i][k] = (M[i][k] - sum) / M[k][k];	//����L����Ĳ���
		}
	}
	
	L = new double *[n];
	for (i = 0; i < n; i++)
		L[i] = new double[n];
	U = new double *[n];
	for (i = 0; i < n; i++)
		U[i] = new double[n];
	L_ = new double *[n];
	for (i = 0; i < n; i++)
		L_[i] = new double[n];
	U_ = new double *[n];
	for (i = 0; i < n; i++)
		U_[i] = new double[n];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			L_[i][j] = U_[i][j] = 0;

	//��L��U���и�ֵ
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i > j) {
				L[i][j] = M[i][j];
				U[i][j] = 0;
			}
			else if (i == j) {
				L[i][j] = 1;
				U[i][j] = M[i][j];
			}
			else {
				L[i][j] = 0;
				U[i][j] = M[i][j];
			}
		}
	}

	for (i = 0; i < n; i++) 
		t *= U[i][i];
	
	if (t == 0.0) {
		cout << "����󲻴���" << endl;
	}
	//�ֱ��������
	if (t) {
		for (i = 0; i < n; i++) {
			U_[i][i] = 1 / U[i][i];
			for (k = i - 1; k >= 0; k--) {
				double sum = 0;
				for (j = k + 1; j <= i; j++)
					sum += U[k][j] * U_[j][i];
				U_[k][i] = -sum / U[k][k];
			}
		}

		for (i = 0; i < n; i++) {
			L_[i][i] = 1;
			for (k = i + 1; k < n; k++) {
				for (j = i; j <= k - 1; j++)
					L_[k][i] -= L[k][j] * L_[j][i];
			}
		}

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				M[i][j] = 0;
				for (k = 0; k < n; k++) {
					M[i][j] += U_[i][k] * L_[k][j];
				}
			}
		}
	}

}