#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<omp.h>
#include<immintrin.h>

int main()
{
    int N = 1e9;
    long long int interval, i;
    double rand_x, rand_y, origin_dist, pi;
    long long int circle_points = 0, square_points = 0;

    // start timer
    clock_t start = clock(), diff;

    square_points = N;
    long long int i = 0, j = 0;
    long long int circle_points = 0;
    double rand_x, rand_y, origin_dist;
    __m256d v_rand_x, v_rand_y, v_origin_dist;
    __m256d v_one = _mm256_set1_pd(1.0);
    long long int v_circle_points = 0;
    for (i = 0; i < N; i += 4)
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

    printf("Estimated value of pi is %f", pi);

    // stop timer
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

    pi = (double)(4 * circle_points) / square_points;
    printf("Estimated value of pi is %f", pi);

    return 0;

}