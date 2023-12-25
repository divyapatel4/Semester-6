#include<stdio.h>
#include<math.h>
#include<stdlib.h>
void transpose_cache(double **src, double **dst, int n)
{
    // Define a block size for cache optimization
    int blocksize = 32;
    int i, j, k, l;

    // Iterate over the matrix in blocks
    for (i = 0; i < n; i += blocksize)
    {
        for (j = 0; j < n; j += blocksize)
        {
            // Transpose the current block
            for (k = i; k < i + blocksize && k < n; ++k)
            {
                for (l = j; l < j + blocksize && l < n; ++l)
                {
                    dst[k][l] = src[l][k];
                }
            }
        }
    }
}

int main(){
    int n = 10;
    double **src = (double **)malloc(n * sizeof(double *));
    double **dst = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        src[i] = (double *)malloc(n * sizeof(double));
        dst[i] = (double *)malloc(n * sizeof(double));
    }
    // random initialization
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            src[i][j] = rand() % 100;
        }
    }

    // print the matrix
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", src[i][j]);
        }
        printf("\n");
    }

    transpose_cache(src, dst, n);

    printf("\n");
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", dst[i][j]);
        }
        printf("\n");
    }

    return 0;
}