
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <immintrin.h>
#define dtype double

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define CLK CLOCK_MONOTONIC

void print_matrix(double **a, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n + 1; j++)
			printf("%lf ", a[i][j]);
		printf("\n");
	}
}

void gaussian_elimination(double **a, int n)
{
	int i, j, k;
	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			double ratio = a[j][i] / a[i][i];
			for (k = i; k < n + 1; k++)
				a[j][k] -= ratio * a[i][k];
		}
	}
}

void back_substitution(double **a, double *x, int n)
{
	int i, j;
	for (i = n - 1; i >= 0; i--)
	{
		x[i] = a[i][n];
		double temp = 0;
		for (j = i + 1; j < n; j++)
			temp += a[i][j] * x[j];
		x[i] -= temp;
		x[i] /= a[i][i];
	}
}

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

	long long int N = atoi(argv[1]); /* size of input array */
	long long int P = atoi(argv[2]); /* number of processors*/
	char *problem_name = "matrix_multiplication";
	char *approach_name = "block";
	FILE *outputFile;

	char outputFileName[50];
	sprintf(outputFileName, "output/%s_%s_%s_%s_output.txt", problem_name, approach_name, argv[1], argv[2]);

	//***************
	long long int i = 0, j = 0;

	// Initializing rand() with seed
	srand(time(NULL));
	//***************
	clock_gettime(CLK, &start_alg); /* Start the algo timer */
	double **A = (double **)malloc(N * sizeof(double *));
	for (i = 0; i < N; i++)
	{
		A[i] = (double *)malloc((N + 1) * sizeof(double));
	}

	double *x = (double *)malloc(N * sizeof(double));

	// Initialize A and x
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N + 1; j++)
		{
			A[i][j] = rand() % 10;
		}
		x[i] = rand() % 10;
	}
	/*----------------------Core algorithm starts here----------------------------------------------*/

	gaussian_elimination(A, N);
	back_substitution(A, x, N);

	print_matrix(A, N);

	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg); /* End the algo timer */

	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;
}
