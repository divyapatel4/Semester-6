#include <stdio.h>
#include <stdlib.h>

void print_matrix(double **a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
            printf("%lf ", a[i][j]);
        printf("\n");
    }
}

void gaussian_elimination(double **a, int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double ratio = a[j][i] / a[i][i];
            for (int k = i; k < n + 1; k++)
                a[j][k] -= ratio * a[i][k];
        }
    }
}

void back_substitution(double **a, double *x, int n)
{
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }
}

int main()
{
    int n;
    printf("Enter the value of n: ");
    scanf("%d", &n);

    double **a = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        a[i] = malloc((n + 1) * sizeof(double));

    double *x = malloc(n * sizeof(double));

    printf("Enter the augmented matrix:\n");
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n + 1; j++)
            scanf("%lf", &a[i][j]);

    printf("Original matrix:\n");
    print_matrix(a, n);

    gaussian_elimination(a, n);

    printf("Matrix after Gaussian elimination:\n");
    print_matrix(a, n);

    back_substitution(a, x, n);

    printf("Solution:\n");
    for (int i = 0; i < n; i++)
        printf("x%d = %lf\n", i + 1, x[i]);

    for (int i = 0; i < n; i++)
        free(a[i]);
    free(a);
    free(x);

    return 0;
}