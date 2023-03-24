
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#define dtype double

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define CLK CLOCK_MONOTONIC

double flatA[2000000];
double flatB[2000000];

struct timespec diff(struct timespec start, struct timespec end);

void convert(double **A, double **B, int n)
{
	int i, j, k;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			flatA[i * n + j] = A[i][j];
			flatB[j * n + i] = B[i][j];
		}
	}
}

void block_matmul(int n, double **A, double **B, double **C)
{
	// Refer optimizedParallelMultiply and make code in format of block_matmul
	int i, j, k, iOff, jOff; // iOff and jOff are used to store the offset of the current block
	dtype tot;				 // tot is used to store the sum of the multiplication of the elements of the current block

	// Convert the matrices A and B to flatA and flatB
	convert(A, B, n);

	for (i = 0; i < n; i++)
	{
		iOff = i * n;
		for (j = 0; j < n; j++)
		{
			jOff = j * n;
			tot = 0;
			for (k = 0; k < n; k++)
			{
				tot += flatA[iOff + k] * flatB[jOff + k];
			}
			C[i][j] = tot;
		}
	}
}

struct timespec diff(struct timespec start, struct timespec end)
{
	struct timespec temp;
	if ((end.tv_nsec - start.tv_nsec) < 0)
	{
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
	}
	else
	{
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}
	return temp;
}

int main(int argc, char *argv[])
{
	struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
	clock_gettime(CLK, &start_e2e);

	/* Check if enough command-line arguments are taken in. */
	if (argc < 3)
	{
		printf("Usage: %s n p \n", argv[0]);
		return -1;
	}

	int N = atoi(argv[1]); /* size of input array */
	int P = atoi(argv[2]); /* number of processors*/
	char *problem_name = "matrix_multiplication";
	char *approach_name = "block";
	FILE *outputFile;

	char outputFileName[50];
	sprintf(outputFileName, "output/%s_%s_%s_%s_output.txt", problem_name, approach_name, argv[1], argv[2]);

	//***************
	// Matrix A of size N*N and Matrix B of size N*N are multiplied to get Matrix C of size N*N

	int i, j, k;

	// using dynamic memory allocation as memory required is very large for large values of N, and stack memory is limited in size so it will cause stack overflow and we will get segmentation fault
	// heap memory is very large in size so we can use it to allocate memory for large values of N

	dtype **A = (dtype **)malloc(N * sizeof(dtype *));
	dtype **B = (dtype **)malloc(N * sizeof(dtype *));
	dtype **C = (dtype **)malloc(N * sizeof(dtype *));

	for (i = 0; i < N; i++)
	{
		A[i] = (dtype *)malloc(N * sizeof(dtype));
		B[i] = (dtype *)malloc(N * sizeof(dtype));
		C[i] = (dtype *)malloc(N * sizeof(dtype));
	}

	

	// Initialize the matrices
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			A[i][j] = 1;
			B[i][j] = 2;
			C[i][j] = 0;
		}
	}

	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/

	block_matmul(N, A, B, C);

	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg); /* End the algo timer */

	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;
}
