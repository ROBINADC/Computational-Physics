#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define N 5	//�����������ά�ȣ���ʱû��ʱ�����ɿɵ�ά��

int main(int argc, char *argv[]) {
	int i, j, k, rank, size, totaltag;
	double A[N][N] = { {0.} }, B[N] = { 0. }, C[N] = { 0. }, D[N] = {0.};
	int major_rank=0; //���������
	double buf[N] = {0.}, starttime, endtime, total = 0.,sum=0.;

	MPI_Status status; 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == major_rank) {   
		//��������A����B
		for (i = 0; i < N; i++) {
			B[i] = rand() / (RAND_MAX + 1.);
			for (j = 0; j < N; j++) {
				A[i][j] = rand() / (RAND_MAX + 1.); 
			}
		}

		//���A,B
		
		printf("Matrix is:\n");
		for(i=0;i<N;i++){
			for (j = 0; j < N; j++) 
				printf("%.3f\t", A[i][j]);
			printf("\n");
		}
		printf("\nVector is:\n");
		for (i = 0; i < N; i++)
			printf("%.3f\n", B[i]);

		//��A���� ɢ������������
		for (i = 1; i < size; i++) {	   
			for (k = i - 1; k < N; k += size - 1) {  
				for (j = 0; j < N; j++) 
					buf[j] = A[k][j]; 
				MPI_Send(buf, N, MPI_DOUBLE, i, k, MPI_COMM_WORLD);
			}
		}
	}

	//��B�㲥����������
	MPI_Bcast(B, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	starttime = MPI_Wtime();	//��¼��ʼʱ��

	if(rank!=major_rank) {	
		for (k = rank-1; k < N; k += size - 1) {	
			//�������̽��ղ����㣬Ȼ���ͻ�������
			MPI_Recv(buf, N, MPI_DOUBLE, major_rank, k, MPI_COMM_WORLD, &status);
			for (j = 0; j < N; j++) 
				total += buf[j] * B[j];
			MPI_Send(&total, 1, MPI_DOUBLE, major_rank, k, MPI_COMM_WORLD);
		}
	}

	if (rank == major_rank) {
		//�����̽���	
		for (i = 0; i < N; i++){
			MPI_Recv(&total, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			totaltag = status.MPI_TAG;
			C[totaltag] = total;
		}
		
		//������
		endtime = MPI_Wtime();
		printf("\nUse %.8fs to complete the calculation.\n", endtime - starttime);
		printf("\nThe result of multiplying is:\n");
		for (i = 0; i < N; i++) {
			printf("%.3f\n", C[i]);
		}

		//��һ��
		for (i = 0; i < N; i++) 
			sum += C[i] * C[i];
		for (i = 0; i < N; i++)
			D[i] = C[i] / sqrt(sum);
		printf("\nNormalized vector:\n");
		for (i = 0; i < N; i++)
			printf("%.3f\n", D[i]);
	}

	MPI_Finalize();
}
