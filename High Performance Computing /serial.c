
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#define dype int
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization("unroll-loops")

#define min(x, y) (((x) < (y)) ? (x) : (y))
//  Using the MONOTONIC clock
#define CLK CLOCK_MONOTONIC

/* Function to compute the difference between two points in time */
struct timespec diff(struct timespec start, struct timespec end);


// Function to multiply two matrices using Strassen's algorithm
void strassen(dype n, dype **A, dype **B, dype **C)
{
    // Base case: if the matrix is 1x1, simply multiply the elements
    if (n == 1)
    {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }

    // Divide the matrices into four submatrices
    dype m = n / 2;

    int i,j;
    // allocate all the matrices dynamically
    dype **A11 = (dype **)malloc(m * sizeof(dype *));
    dype **A12 = (dype **)malloc(m * sizeof(dype *));
    dype **A21 = (dype **)malloc(m * sizeof(dype *));
    dype **A22 = (dype **)malloc(m * sizeof(dype *));
    dype **B11 = (dype **)malloc(m * sizeof(dype *));
    dype **B12 = (dype **)malloc(m * sizeof(dype *));
    dype **B21 = (dype **)malloc(m * sizeof(dype *));
    dype **B22 = (dype **)malloc(m * sizeof(dype *));
    dype **C11 = (dype **)malloc(m * sizeof(dype *));
    dype **C12 = (dype **)malloc(m * sizeof(dype *));
    dype **C21 = (dype **)malloc(m * sizeof(dype *));
    dype **C22 = (dype **)malloc(m * sizeof(dype *));
    dype **P1 = (dype **)malloc(m * sizeof(dype *));
    dype **P2 = (dype **)malloc(m * sizeof(dype *));
    dype **P3 = (dype **)malloc(m * sizeof(dype *));
    dype **P4 = (dype **)malloc(m * sizeof(dype *));
    dype **P5 = (dype **)malloc(m * sizeof(dype *));
    dype **P6 = (dype **)malloc(m * sizeof(dype *));
    dype **P7 = (dype **)malloc(m * sizeof(dype *));

    for ( i = 0; i < m; i++)
    {
        A11[i] = (dype *)malloc(m * sizeof(dype));
        A12[i] = (dype *)malloc(m * sizeof(dype));
        A21[i] = (dype *)malloc(m * sizeof(dype));
        A22[i] = (dype *)malloc(m * sizeof(dype));
        B11[i] = (dype *)malloc(m * sizeof(dype));
        B12[i] = (dype *)malloc(m * sizeof(dype));
        B21[i] = (dype *)malloc(m * sizeof(dype));
        B22[i] = (dype *)malloc(m * sizeof(dype));
        C11[i] = (dype *)malloc(m * sizeof(dype));
        C12[i] = (dype *)malloc(m * sizeof(dype));
        C21[i] = (dype *)malloc(m * sizeof(dype));
        C22[i] = (dype *)malloc(m * sizeof(dype));
        P1[i] = (dype *)malloc(m * sizeof(dype));
        P2[i] = (dype *)malloc(m * sizeof(dype));
        P3[i] = (dype *)malloc(m * sizeof(dype));
        P4[i] = (dype *)malloc(m * sizeof(dype));
        P5[i] = (dype *)malloc(m * sizeof(dype));
        P6[i] = (dype *)malloc(m * sizeof(dype));
        P7[i] = (dype *)malloc(m * sizeof(dype));
    }

    for ( i = 0; i < m; i++)
    {
        for ( j = 0; j < m; j++)
        {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + m];
            A21[i][j] = A[i + m][j];
            A22[i][j] = A[i + m][j + m];
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + m];
            B21[i][j] = B[i + m][j];
            B22[i][j] = B[i + m][j + m];
        }
    }

    // Calculate seven products recursively
    strassen(m, A11, B11, P1);
    strassen(m, A12, B21, P2);
    strassen(m, A11, B12, P3);
    strassen(m, A12, B22, P4);
    strassen(m, A21, B11, P5);
    strassen(m, A22, B21, P6);
    strassen(m, A21, B12, P7);
    strassen(m, A22, B22, C22);

    // Calculate the submatrices of the resulting matrix
    for ( i = 0; i < m; i++)
    {
        for ( j = 0; j < m; j++)
        {
            C11[i][j] = P1[i][j] + P2[i][j];
            C12[i][j] = P3[i][j] + P4[i][j];
            C21[i][j] = P5[i][j] + P6[i][j];
            C22[i][j] = P7[i][j] + C22[i][j];
        }
    }

    // Combine the submatrices into the final matrix
    for ( i = 0; i < m; i++)
    {
        for ( j = 0; j < m; j++)
        {
            C[i][j] = C11[i][j];
            C[i][j + m] = C12[i][j];
            C[i + m][j] = C21[i][j];
            C[i + m][j + m] = C22[i][j];
        }
    }
    // free all the allocated memory
    for ( i = 0; i < m; i++)
    {
        free(A11[i]);
        free(A12[i]);
        free(A21[i]);
        free(A22[i]);
        free(B11[i]);
        free(B12[i]);
        free(B21[i]);
        free(B22[i]);
        free(C11[i]);
        free(C12[i]);
        free(C21[i]);
        free(C22[i]);
        free(P1[i]);
        free(P2[i]);
        free(P3[i]);
        free(P4[i]);
        free(P5[i]);
        free(P6[i]);
        free(P7[i]);
    }
    free(A11); free(A12); free(A21); free(A22); free(B11); free(B12); free(B21); free(B22); free(C11); free(C12); free(C21); free(C22); free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7);

            
}

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
    // Problem - 1 Matrix Multiplication

    // Matrix A of size N*N and Matrix B of size N*N are multiplied to get Matrix C of size N*N

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

    // Matrix multiplication using Strassen's algorithm
    strassen(N, A, B, C);

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
