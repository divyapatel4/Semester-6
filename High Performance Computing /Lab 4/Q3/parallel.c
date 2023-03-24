
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))
//  Using the MONOTONIC clock
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
	/* Should start before anything else */
	clock_gettime(CLK, &start_e2e);

	/* Check if enough command-line arguments are taken in. */
	if (argc < 3)
	{
		printf("Usage: %s n p \n", argv[0]);
		return -1;
	}

	int N = atoi(argv[1]); /* size of input array */
	int P = atoi(argv[2]); /* number of threads */
	char *problem_name = "matrix_multiplication";
	char *approach_name = "block";
	//	char buffer[10];
	//	FILE* inputFile;
	FILE *outputFile;
	//	inputFile = fopen(argv[3],"r");

	char outputFileName[50];
	sprintf(outputFileName, "output/%s_%s_%s_%s_output.txt", problem_name, approach_name, argv[1], argv[2]);

	//***************
	// Problem -3
	double *arr = malloc(N * sizeof(double));
	int i = 0;
	omp_set_num_threads(P);

	#pragma omp for
	for (i = 0; i < N; i++)
	{
		arr[i] = N - i;
	}
	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/

	void quick_sort_parallel(double *arr, int left, int right)
	{
		int i = left;
		int j = right;
		double tmp;
		double pivot = arr[(left + right) / 2];

		/* partition */
		while (i <= j)
		{
			while (arr[i] < pivot)
				i++;
			while (arr[j] > pivot)
				j--;
			if (i <= j)
			{
				tmp = arr[i];
				arr[i] = arr[j];
				arr[j] = tmp;
				i++;
				j--;
			}
		}

		/* recursion */
		#pragma omp task firstprivate(arr, left, j)
		{
			if (left < j)
			{
				quick_sort_parallel(arr, left, j);
			}
		}
		#pragma omp task firstprivate(arr, right, i)
		{
			if (i < right)
			{
				quick_sort_parallel(arr, i, right);
			}
		}
	}

	void quick_sort(double *arr, int N, int P)
	{
		omp_set_num_threads(P);
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				quick_sort_parallel(arr, 0, N - 1);
			}
		}
	}

	quick_sort(arr, N, P);

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
