
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include<immintrin.h>
#define dtype double

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

	long long int N = atoi(argv[1]); /* size of input array */
	long long int P = atoi(argv[2]); /* number of processors*/
	char *problem_name = "matrix_multiplication";
	char *approach_name = "block";
	FILE *outputFile;

	char outputFileName[50];
	sprintf(outputFileName, "output/%s_%s_%s_%s_output.txt", problem_name, approach_name, argv[1], argv[2]);

	//***************
	long long int interval, i = 0, j = 0;
	// double rand_x, rand_y, origin_dist, pi;
	double pi;
	unsigned int seed = 42;
	long long int circle_points = 0, square_points = 0;

	// Initializing rand() with seed
	srand(time(NULL));

	//***************

	clock_gettime(CLK, &start_alg); /* Start the algo timer */
	long long INTERVAL = N;

	/*----------------------Core algorithm starts here----------------------------------------------*/

	double rand_x, rand_y, origin_dist;
	__m256d v_rand_x, v_rand_y, v_origin_dist;
	__m256d v_one = _mm256_set1_pd(1.0);
	long long int v_circle_points = 0;
	for ( i = 0; i < N; i += 4)
	{
		// Randomly generated x and y values
		v_rand_x = _mm256_set_pd((double)rand_r(&seed) / RAND_MAX,
								 (double)rand_r(&seed) / RAND_MAX,
								 (double)rand_r(&seed) / RAND_MAX,
								 (double)rand_r(&seed) / RAND_MAX);
		v_rand_y = _mm256_set_pd((double)rand_r(&seed) / RAND_MAX,
								 (double)rand_r(&seed) / RAND_MAX,
								 (double)rand_r(&seed) / RAND_MAX,
								 (double)rand_r(&seed) / RAND_MAX);
		// Distance between (x, y) from the origin
		v_origin_dist = _mm256_add_pd(_mm256_mul_pd(v_rand_x, v_rand_x), _mm256_mul_pd(v_rand_y, v_rand_y));
		// Check if the point is inside the circle
		v_circle_points += _mm_popcnt_u32(_mm256_movemask_pd(_mm256_cmp_pd(v_origin_dist, v_one, _CMP_LE_OS)));
	}
	circle_points += v_circle_points;
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
