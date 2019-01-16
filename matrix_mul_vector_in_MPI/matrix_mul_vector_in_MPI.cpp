#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define N 5	//方阵和向量的维度，暂时没有时间做成可调维度

int main(int argc, char *argv[]) {
	int i, j, k, rank, size, totaltag;
	double A[N][N] = { {0.} }, B[N] = { 0. }, C[N] = { 0. }, D[N] = {0.};
	int major_rank=0; //主进程序号
	double buf[N] = {0.}, starttime, endtime, total = 0.,sum=0.;

	MPI_Status status; 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == major_rank) {   
		//创建矩阵A向量B
		for (i = 0; i < N; i++) {
			B[i] = rand() / (RAND_MAX + 1.);
			for (j = 0; j < N; j++) {
				A[i][j] = rand() / (RAND_MAX + 1.); 
			}
		}

		//输出A,B
		
		printf("Matrix is:\n");
		for(i=0;i<N;i++){
			for (j = 0; j < N; j++) 
				printf("%.3f\t", A[i][j]);
			printf("\n");
		}
		printf("\nVector is:\n");
		for (i = 0; i < N; i++)
			printf("%.3f\n", B[i]);

		//将A按行 散播到其他进程
		for (i = 1; i < size; i++) {	   
			for (k = i - 1; k < N; k += size - 1) {  
				for (j = 0; j < N; j++) 
					buf[j] = A[k][j]; 
				MPI_Send(buf, N, MPI_DOUBLE, i, k, MPI_COMM_WORLD);
			}
		}
	}

	//将B广播到其他进程
	MPI_Bcast(B, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	starttime = MPI_Wtime();	//记录开始时间

	if(rank!=major_rank) {	
		for (k = rank-1; k < N; k += size - 1) {	
			//其他进程接收并计算，然后发送回主进程
			MPI_Recv(buf, N, MPI_DOUBLE, major_rank, k, MPI_COMM_WORLD, &status);
			for (j = 0; j < N; j++) 
				total += buf[j] * B[j];
			MPI_Send(&total, 1, MPI_DOUBLE, major_rank, k, MPI_COMM_WORLD);
		}
	}

	if (rank == major_rank) {
		//主进程接收	
		for (i = 0; i < N; i++){
			MPI_Recv(&total, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			totaltag = status.MPI_TAG;
			C[totaltag] = total;
		}
		
		//输出结果
		endtime = MPI_Wtime();
		printf("\nUse %.8fs to complete the calculation.\n", endtime - starttime);
		printf("\nThe result of multiplying is:\n");
		for (i = 0; i < N; i++) {
			printf("%.3f\n", C[i]);
		}

		//归一化
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
