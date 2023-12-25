
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#define dype double

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define CLK CLOCK_MONOTONIC

struct timespec diff(struct timespec start, struct timespec end);

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
	// Problem - 1 Matrix Multiplication
	int i, j, k;
	// using dynamic memory allocation as memory required is very large for large values of N, and stack memory is limited in size so it will cause stack overflow and we will get segmentation fault
	// heap memory is very large in size so we can use it to allocate memory for large values of N
	dype **A = (dype **)malloc(N * sizeof(dype *));
	dype **B = (dype **)malloc(N * sizeof(dype *));
	dype **C = (dype **)malloc(N * sizeof(dype *));

	for (i = 0; i < N; i++)
	{
		A[i] = (dype *)malloc(N * sizeof(dype));
		B[i] = (dype *)malloc(N * sizeof(dype));
		C[i] = (dype *)malloc(N * sizeof(dype));
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
	dype **B_transpose = (dype **)malloc(N * sizeof(dype *));
	for (i = 0; i < N; i++)
	{
		B_transpose[i] = (dype *)malloc(N * sizeof(dype));
	}

	// Calculate B_transpose
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			B_transpose[i][j] = B[j][i];
		}
	}

	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/

	// Matrix Multiplication by calculating B_transpose and then multiplying it with A

	// Multiply A and B_transpose
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < N; k++)
			{
				C[i][j] += A[i][k] * B_transpose[j][k];
			}
		}
	}

	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg); /* End the algo timer */

	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;
}
