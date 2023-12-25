#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<omp.h>
// #include<immintrin.h>

int main()
{
    long long int N = 1e10;
    long long int interval, i,j;
    double rand_x, rand_y, origin_dist, pi;
    long long int circle_points = 0, square_points = 0;

    // start timer
    clock_t start = clock(), diff;
    unsigned int seed = 0;
    square_points = N;
	#pragma omp parallel
	{
		unsigned int seed = (unsigned)time(NULL)*omp_get_thread_num();
		#pragma omp for private(i, rand_x, rand_y, origin_dist) schedule(static) reduction(+:circle_points) 
		for (i = 0; i < N; i++)
		{
			// Randomly generated x and y values
			rand_x = (rand_r(&seed)/(double)RAND_MAX);
			rand_y = (rand_r(&seed)/(double)RAND_MAX);

			origin_dist = rand_x * rand_x + rand_y * rand_y;

			if (origin_dist <= 1)
				circle_points++;
		}
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