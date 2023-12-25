#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
// #include<omp.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))

void transpose(double **src, double **dst, int n)
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

int main()
{
    int n = 8;
    double **a = (double **)malloc(n * sizeof(double *));
    double **q = (double **)malloc(n * sizeof(double *));
    double **r = (double **)malloc(n * sizeof(double *));
    
    for (int i = 0; i < n; i++)
    {
        a[i] = (double *)malloc(n * sizeof(double));
        q[i] = (double *)malloc(n * sizeof(double));
        r[i] = (double *)malloc(n * sizeof(double));
    }
    a[0][0] = 6.0; a[0][1] = 4.0; a[0][2] = 5.0; a[0][3] = 3.0; a[0][4] = 3.0; a[0][5] = 2.0; a[0][6] = 9.0; a[0][7] = 1.0;
    a[1][0] = 0.0; a[1][1] = 0.0; a[1][2] = 0.0; a[1][3] = 2.0; a[1][4] = 5.0; a[1][5] = 8.0; a[1][6] = 4.0; a[1][7] = 1.0;
    a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 4.0; a[2][3] = 6.0; a[2][4] = 7.0; a[2][5] = 8.0; a[2][6] = 5.0; a[2][7] = 2.0;
    a[3][0] = 0.0; a[3][1] = 0.0; a[3][2] = 0.0; a[3][3] = 0.0; a[3][4] = 0.0; a[3][5] = 0.0; a[3][6] = 9.0; a[3][7] = 6.0;
    a[4][0] = 0.0; a[4][1] = 0.0; a[4][2] = 0.0; a[4][3] = 0.0; a[4][4] = 1.0; a[4][5] = 2.0; a[4][6] = 8.0; a[4][7] = 6.0;
    a[5][0] = 0.0; a[5][1] = 4.0; a[5][2] = 6.0; a[5][3] = 7.0; a[5][4] = 8.0; a[5][5] = 4.0; a[5][6] = 4.0; a[5][7] = 6.0;
    a[6][0] = 0.0; a[6][1] = 0.0; a[6][2] = 0.0; a[6][3] = 0.0; a[6][4] = 0.0; a[6][5] = 0.0; a[6][6] = 0.0; a[6][7] = 5.0;
    a[7][0] = 0.0; a[7][1] = 0.0; a[7][2] = 0.0; a[7][3] = 0.0; a[7][4] = 0.0; a[7][5] = 4.0; a[7][6] = 8.0; a[7][7] = 3.0;

    int N = n,i,j,k;
    double s;

    double **a_transpose = (double **)malloc(N * sizeof(double *)), **q_transpose = (double **)malloc(N * sizeof(double *));
    for (int row = 0; row < N; row++)
    {
        a_transpose[row] = (double *)malloc(N * sizeof(double));
        q_transpose[row] = (double *)malloc(N * sizeof(double));
    }

    transpose(a, a_transpose, N);

    // #pragma omp parallel for private(i, j, k)
    for (i = 0; i < N; i++)
    {
        // Compute the diagonal element of R
        double sum = 0;
        #pragma omp simd reduction(+: sum)
        for (j = 0; j < N; j++)
        {
            sum += a_transpose[i][j] * a_transpose[i][j];
        }
        r[i][i] = sqrt(sum);

        // Compute the i-th column of Q
        #pragma omp simd
        for (j = 0; j < N; j++)
        {
            q_transpose[i][j] = a_transpose[i][j] / r[i][i];
        }

        // Update the remaining columns of R and A_T
        #pragma omp parallel for private(j, k)
        for (j = i + 1; j < N; j++)
        {
            // Compute the off-diagonal elements of R
            double sum = 0;
            #pragma omp simd reduction(+: sum)
            for (k = 0; k < N; k++)
            {
                sum += q_transpose[i][k] * a_transpose[j][k];
            }
            r[i][j] = sum;

            // Update A_T by subtracting the projection of column j onto column i
            #pragma omp simd
            for (k = 0; k < N; k++)
            {
                a_transpose[j][k] -= sum * q_transpose[i][k];
            }
        }
    }

    // Transpose the Q matrix to obtain an orthogonal matrix
    transpose(q_transpose, q, N);

    // print the result
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", q[i][j]);
        }
        printf("\n");
    }
    printf("\n");   

    // Print R
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", r[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < n; i++)
    {
        free(a[i]);
        free(q[i]);
        free(r[i]);
    }
    free(a);
    free(q);
    free(r);
}
