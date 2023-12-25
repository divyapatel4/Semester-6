#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <immintrin.h>
#define dtype double
#define SEED 35791246


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

const int nThreads = 16; // number of threads to use
#define CACHE_LINE_SIZE 64
#define PADDING (CACHE_LINE_SIZE / sizeof(unsigned int))

unsigned int seeds[16 * PADDING];

void seedThreads()
{
	int i;
	for ( i = 0; i < nThreads; i++)
	{
		unsigned int seed = (unsigned)time(NULL);
		seeds[i * PADDING] = (seed & 0xFFFFFFF0) | (i + 1);
	}
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
	long long int interval, i;
	double rand_x, rand_y, origin_dist, pi;
	long long int circle_points = 0, square_points = 0;

	// Initializing rand() with seed
	srand(time(NULL));

	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */
	long long INTERVAL = N;

	omp_set_num_threads(P);
	seedThreads();
	/*----------------------Core algorithm starts here----------------------------------------------*/
    square_points = N;
	
	#pragma omp parallel
	{
		int tid = omp_get_thread_num(); 
		unsigned int seed = seeds[tid * PADDING]; 
		srand(seed); 
		double rand_x, rand_y, origin_dist; 
		#pragma omp for private(i) reduction(+ : circle_points) schedule(static) 
		for (i = 0; i < N; i++) { 
			// Randomly generated x and y values 
			rand_x = (double)rand_r(&seed) / RAND_MAX; 
			rand_y = (double)rand_r(&seed) / RAND_MAX; 
			
			// Distance between (x, y) from the origin 
			origin_dist = rand_x * rand_x + rand_y * rand_y; 
			if (origin_dist <= 1) 
				circle_points++; 
		} 
	}

	pi = (double)(4 * circle_points) / square_points;

	// printf("Estimated value of pi is : %f \n", pi);
	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg); /* End the algo timer */

	clock_gettime(CLK, &end_e2e);

	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;
}

