#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<time.h>
#include<string.h>
#include<stdlib.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))
//  Using the MONOTONIC clock 
#define CLK CLOCK_MONOTONIC

/* Function to compute the difference between two points in time */
struct timespec diff(struct timespec start, struct timespec end);

struct timespec diff(struct timespec start, struct timespec end){
	struct timespec temp;
	if((end.tv_nsec-start.tv_nsec)<0){
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

double l2_norm(double *v, int n)
{
	double s = 0.0;
	// #pragma omp parallel for reduction(+:s)
	for (int i = 0; i < n; i++)
	{
		s += v[i] * v[i];
	}
	return sqrt(s);
}

void modifiedGS(double *A, int m, int n, double *Q, double *R)
{
	// Allocate memory for temporary matrix
	double *V = (double *)malloc(sizeof(double) * m * n);
	memcpy(V, A, sizeof(double) * m * n);

	// Set R to zero
	memset(R, 0, sizeof(double) * n * n);

	// Iterate over columns
	for (int k = 0; k < n; k++)
	{
		// Compute R_kk and Q_k
		double norm = 0.0;
		// #pragma omp parallel for reduction(+ \
                                   : norm)
		for (int i = 0; i < m; i++)
		{
			norm += V[i * n + k] * V[i * n + k];
		}
		norm = sqrt(norm);
		R[k * n + k] = norm;
		// #pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			Q[i * n + k] = V[i * n + k] / norm;
		}

		// Update V and R
		for (int j = k + 1; j < n; j++)
		{
			double dot = 0.0;
			// #pragma omp parallel for reduction(+ \
                                   : dot)
			for (int i = 0; i < m; i++)
			{
				dot += Q[i * n + k] * V[i * n + j];
			}
			R[k * n + j] = dot;
			// #pragma omp parallel for
			for (int i = 0; i < m; i++)
			{
				V[i * n + j] -= dot * Q[i * n + k];
			}
		}
	}

	// Free memory
	free(V);
}

void blockQR(double *A, int m, int n, int b, double *Q, double *R)
{
	// Allocate memory for temporary matrices
	double *V = (double *)malloc(sizeof(double) * m * b);
	double *W = (double *)malloc(sizeof(double) * n * b);

	// Set R to zero
	memset(R, 0, sizeof(double) * n * n);

	// Iterate over blocks
	for (int k = 0; k < n; k += b)
	{
		int bs = min(b, n - k);

		// Copy block from A to V
		// #pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			memcpy(V + i * bs, A + i * n + k, sizeof(double) * bs);
		}

		// Update V and R
		for (int j = 0; j < k; j += b)
		{
			int bs2 = min(b, n - j);

			// Compute W = Q_j^T * V
			// #pragma omp parallel for
			for (int p = 0; p < bs2; p++)
			{
				for (int q = 0; q < bs; q++)
				{
					double sum = 0.0;
					for (int i = 0; i < m; i++)
					{
						sum += Q[i * n + j + p] * V[i * bs + q];
					}
					W[p * bs + q] = sum;
				}
			}

			// Update R
			// #pragma omp parallel for
			for (int p = 0; p < bs2; p++)
			{
				for (int q = 0; q < bs; q++)
				{
					R[(j + p) * n + k + q] = W[p * bs + q];
				}
			}

			// Update V
			// #pragma omp parallel for
			for (int i = 0; i < m; i++)
			{
				for (int q = 0; q < bs; q++)
				{
					double sum = 0.0;
					for (int p = 0; p < bs2; p++)
					{
						sum += Q[i * n + j + p] * W[p * bs + q];
					}
					V[i * bs + q] -= sum;
				}
			}
		}

		// Compute QR decomposition of V
		double *Q2 = (double *)malloc(sizeof(double) * m * bs);
		double *R2 = (double *)malloc(sizeof(double) * bs * bs);
		modifiedGS(V, m, bs, Q2, R2);

		// Update Q and R
		// #pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			memcpy(Q + i * n + k, Q2 + i * bs, sizeof(double) * bs);
		}
		// #pragma omp parallel for
		for (int p = 0; p < bs; p++)
		{
			memcpy(R + (k + p) * n + k, R2 + p * bs, sizeof(double) * bs);
		}

		// Free memory
		free(Q2);
		free(R2);
	}

	// Free memory
	free(V);
	free(W);
}

void qr_decomposition_new3(double **a, int m, int n, int b, double **q, double **r)
{
	// Convert 2D arrays to 1D arrays
	double *a_1d = (double *)malloc(sizeof(double) * m * n);
	double *q_1d = (double *)malloc(sizeof(double) * m * n);
	double *r_1d = (double *)malloc(sizeof(double) * n * n);

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a_1d[i * n + j] = a[i][j];
		}
	}

	// Call blockQR function
	blockQR(a_1d, m, n, b, q_1d, r_1d);

	// Convert 1D arrays back to 2D arrays
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			q[i][j] = q_1d[i * n + j];
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			r[i][j] = r_1d[i * n + j];
		}
	}

	// Free allocated memory
	free(a_1d);
	free(q_1d);
	free(r_1d);
}

void classicalGS(double *A_current, double *A_T, int P, double *Q_current, double *R_current)
{
    double *v_col = (double *)malloc(sizeof(double) * P);
    int col, row, row_;
    double result;
    for (col = 0; col < P; col++)
    {
        memcpy(v_col, A_T + col * P, sizeof(double) * P);
        for (row = 0; row < col; row++)
        {
            result = 0.0;
            #pragma omp parallel for reduction(+:result)
            for (row_ = 0; row_ < P; row_++)
            {
                result += Q_current[row_ * P + row] * A_T[col * P + row_];
            }
            R_current[row * P + col] = result;
            #pragma omp parallel for
            for (row_ = 0; row_ < P; row_++)
            {
                v_col[row_] -= R_current[row * P + col] * Q_current[row_ * P + row];
            }
        }
        R_current[col * P + col] = l2_norm(v_col, P);
        #pragma omp parallel for
        for (row = 0; row < P; row++)
        {
            Q_current[row * P + col] = v_col[row] / R_current[col * P + col];
        }
    }
    free(v_col);
}

void qr_decomposition_new(double **a, double **q, double **r, int n)
{
	// Allocate memory for the input matrices in the format required by classicalGS
	double *A_current = (double *)malloc(sizeof(double) * n * n);
	double *A_T = (double *)malloc(sizeof(double) * n * n);
	double *Q_current = (double *)malloc(sizeof(double) * n * n);
	double *R_current = (double *)malloc(sizeof(double) * n * n);

	// Copy the input matrix to A_current and A_T
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A_current[i * n + j] = a[i][j];
			A_T[j * n + i] = a[i][j];
		}
	}

	// Call the classicalGS function
	classicalGS(A_current, A_T, n, Q_current, R_current);

	// Copy the results from Q_current and R_current to the q and r matrices
	#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			q[i][j] = Q_current[i * n + j];
			r[i][j] = R_current[i * n + j];
		}
	}

	// Free the memory allocated for the input matrices
	free(A_current);
	free(A_T);
	free(Q_current);
	free(R_current);
}

int main(int argc, char* argv[])
{
	struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
	/* Should start before anything else */
	clock_gettime(CLK, &start_e2e);

	/* Check if enough command-line arguments are taken in. */
	if(argc < 3){
		printf( "Usage: %s N p \n", argv[0] );
		return -1;
	}

	int N=atoi(argv[1]);	/* size of input array */
	int P=atoi(argv[2]);	/* number of processors*/
	char *problem_name = "matrix_multiplication";
	char *approach_name = "block";
//	char buffer[10];
//	FILE* inputFile;
	FILE* outputFile;
	//	inputFile = fopen(argv[3],"r");

	char outputFileName[50];		
	// sprintf(outputFileName,"output/%s_%s_%s_%s_output.txt",problem_name,approach_name,argv[1],argv[2]);
	snprintf(outputFileName, sizeof(outputFileName), "output/%s_%s_%s_%s_output.txt", problem_name, approach_name, argv[1], argv[2]);

	//***************
	
    int i, j, k;
	// int N = 8;
    // double a[8][8];
    double **a = (double **)malloc(N * sizeof(double *));
    for(i = 0; i < N; i++) {
        a[i] = (double *)malloc(N * sizeof(double));
    }
    double **q = (double **)malloc(N * sizeof(double *));
    for(i = 0; i < N; i++) {
        q[i] = (double *)malloc(N * sizeof(double));
    }
    double **r = (double **)malloc(N * sizeof(double *));
    for(i = 0; i < N; i++) {
        r[i] = (double *)malloc(N * sizeof(double));
    }

    srand(time(0));
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
            a[i][j] = rand()%10;
        }
    }

    double s;
	
	
	//***************


	clock_gettime(CLK, &start_alg);	/* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/

	qr_decomposition_new3(a, N, N, 32, q, r);

	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg);	/* End the algo timer */
	/* Ensure that only the algorithm is present between these two
	   timers. Further, the whole algorithm should be present. */


	/* Should end before anything else (printing comes later) */
	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;

}
