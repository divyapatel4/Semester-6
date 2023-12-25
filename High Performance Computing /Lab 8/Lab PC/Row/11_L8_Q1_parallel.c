
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#define dype float
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define CLK CLOCK_MONOTONIC


void rowWiseGaussianEliminationParallel(double **a, int N) {
	int i, j, k;
	for ( k = 0; k < N - 1; k++) {
        #pragma omp parallel for collapse(2)
	
		for ( i = k + 1; i < N; i++) {
			
			for ( j = k + 1; j <= N; j++) {
                a[i][j] -= (a[i][k] / a[k][k]) * a[k][j];
            }
        }
    }
}

void backSubstitutionParallel(double **a, double x[], int N) {
    x[N-1] = a[N-1][N] / a[N-1][N-1];
	int i, j;
    for ( i = N - 2; i >= 0; i--) {
        double sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for ( j = i + 1; j < N; j++) {
            sum += a[i][j] * x[j];
        }
        x[i] = (a[i][N] - sum) / a[i][i];
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

#define INTERVAL 10000
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
	int i, j, k;
	double **a = (double **)malloc(N * sizeof(double *));
	for (i = 0; i < N; i++)
	{
		a[i] = (double *)malloc((N + 1) * sizeof(double));
	}
	double *x = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		x[i] = i + 1;
	}
	for (i = 0; i < N; i++)
	{
		a[i][N] = 0;
		for (j = 0; j < N; j++)
		{
			a[i][j] = 1.0 / (i + j + 1);
			a[i][N] += a[i][j] * x[j];
		}
	}

	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/
	omp_set_num_threads(P);

	rowWiseGaussianEliminationParallel(a, N);
	backSubstitutionParallel(a, x, N);

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
