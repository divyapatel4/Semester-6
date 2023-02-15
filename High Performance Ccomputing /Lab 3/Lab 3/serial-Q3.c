
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
// #pragma GCC optimize("Ofast")
// #pragma GCC target("avx,avx2,fma")
// #pragma GCC optimization("unroll-loops")

#define min(x, y) (((x) < (y)) ? (x) : (y))
//  Using the MONOTONIC clock
#define CLK CLOCK_MONOTONIC

#define dype int

#define BLOCK_SIZE 4

void add(int n, int A[n][n], int B[n][n], int C[n][n])
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			C[i][j] = A[i][j] + B[i][j];
}

void subtract(int n, int A[n][n], int B[n][n], int C[n][n])
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			C[i][j] = A[i][j] - B[i][j];
}

void strassen_morton(int n, int A[n][n], int B[n][n], int C[n][n])
{
	if (n <= BLOCK_SIZE)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				C[i][j] = 0;
				for (int k = 0; k < n; k++)
					C[i][j] += A[i][k] * B[k][j];
			}
		return;
	}

	int half = n / 2;
	int(*A11)[half] = A, (*A12)[half] = A + half, (*A21)[half] = A + n * half, (*A22)[half] = A + n * half + half;
	int(*B11)[half] = B, (*B12)[half] = B + half, (*B21)[half] = B + n * half, (*B22)[half] = B + n * half + half;
	int(*C11)[half] = C, (*C12)[half] = C + half, (*C21)[half] = C + n * half, (*C22)[half] = C + n * half + half;
	int(*P1)[half] = malloc(half * half * sizeof(int)), (*P2)[half] = malloc(half * half * sizeof(int));
	int(*P3)[half] = malloc(half * half * sizeof(int)), (*P4)[half] = malloc(half * half * sizeof(int));
	int(*P5)[half] = malloc(half * half * sizeof(int)), (*P6)[half] = malloc(half * half * sizeof(int));
	int(*P7)[half] = malloc(half * half * sizeof(int));
	int(*T1)[half] = malloc(half * half * sizeof(int)), (*T2)[half] = malloc(half * half * sizeof(int));

	subtract(half, B12, B22, T1);		// T1 = B12 - B22
	strassen_morton(half, A11, T1, P1); // P1 = A11 * T1
	add(half, A11, A12, T1);			// T1 = A11 + A12
	strassen_morton(half, T1, B22, P2); // P2 = T1 * B22

	add(half, A21, A22, T1);			// T1 = A21 + A22
	strassen_morton(half, T1, B11, P3); // P3 = T1 * B11

	subtract(half, B21, B11, T1);		// T1 = B21 - B11
	strassen_morton(half, A22, T1, P4); // P4 = A22 * T1

	add(half, A11, A22, T1);		   // T1 = A11 + A22
	add(half, B11, B22, T2);		   // T2 = B11 + B22
	strassen_morton(half, T1, T2, P5); // P5 = T1 * T2

	subtract(half, A12, A22, T1);	   // T1 = A12 - A22
	add(half, B21, B22, T2);		   // T2 = B21 + B22
	strassen_morton(half, T1, T2, P6); // P6 = T1 * T2

	subtract(half, A11, A21, T1);	   // T1 = A11 - A21
	add(half, B11, B12, T2);		   // T2 = B11 + B12
	strassen_morton(half, T1, T2, P7); // P7 = T1 * T2

	add(half, P5, P4, T1);		// T1 = P5 + P4
	subtract(half, T1, P2, T2); // T2 = T1 - P2
	add(T2, P6, C11);			// C11 = T2 + P6

	add(half, P1, P2, C12); // C12 = P1 + P2

	add(half, P3, P4, C21); // C21 = P3 + P4

	subtract(half, P5, P3, T1); // T1 = P5 - P3
	subtract(half, T1, P7, T2); // T2 = T1 - P7
	add(T2, P6, C22);			// C22 = T2 + P6

	free(P1);
	free(P2);
	free(P3);
	free(P4);
	free(P5);
	free(P6);
	free(P7);
	free(T1);
	free(T2);
}

/* Function to compute the difference between two points in time */
struct timespec diff(struct timespec start, struct timespec end);

/*
   Function to computes the difference between two time instances

   Taken from - http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/

   Further reading:
http://stackoverflow.com/questions/6749621/how-to-create-a-high-resolution-timer-in-linux-to-measure-program-performance
http://stackoverflow.com/questions/3523442/difference-between-clock-realtime-and-clock-monotonic
 */
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
	/* Should start before anything else */
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
	//	char buffer[10];
	//	FILE* inputFile;
	FILE *outputFile;
	//	inputFile = fopen(argv[3],"r");

	char outputFileName[50];
	sprintf(outputFileName, "output/%s_%s_%s_%s_output.txt", problem_name, approach_name, argv[1], argv[2]);

	//***************
	// Problem - 3 Matrix Multiplication using Strassen Algorithm

	int i, j, k;
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

	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/
	// Matrix Multiplication by using block approach with strassen morton algorithm

	int block_size = 4;


	strassen_morton(N, A, B, C);
	
	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg); /* End the algo timer */
	/* Ensure that only the algorithm is present between these two
	   timers. Further, the whole algorithm should be present. */

	/* Should end before anything else (printing comes later) */
	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	/* problem_name,approach_name,n,p,e2e_sec,e2e_nsec,alg_sec,alg_nsec
	   Change problem_name to whatever problem you've been assigned
	   Change approach_name to whatever approach has been assigned
	   p should be 0 for serial codes!!
	 */
	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;
}
