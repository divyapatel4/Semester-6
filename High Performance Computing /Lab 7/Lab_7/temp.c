#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<omp.h>
// #include<immintrin.h>

int main()
{
    int N = 1e9;
    long long int interval, i,j;
    double rand_x, rand_y, origin_dist, pi;
    long long int circle_points = 0, square_points = 0;

    // start timer
    clock_t start = clock(), diff;
    unsigned int seed = 0;
    square_points = N;
	#pragma omp parallel for private(i, rand_x, rand_y, origin_dist) reduction(+ : circle_points) schedule(static)
	for (i = 0; i < N; i++)
	{
		// Randomly generated x and y values
		rand_x = (double)rand_r(&seed) / RAND_MAX;
		rand_y = (double)rand_r(&seed) / RAND_MAX;

		// Distance between (x, y) from the origin
		origin_dist = rand_x * rand_x + rand_y * rand_y;

		if (origin_dist <= 1)
			circle_points++;
	}

    pi = (double)(4 * circle_points) / square_points;

    printf("Estimated value of pi is %f", pi);

    // stop timer
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

    pi = (double)(4 * circle_points) / square_points;
    printf("Estimated value of pi is %f", pi);

    return 0;

}