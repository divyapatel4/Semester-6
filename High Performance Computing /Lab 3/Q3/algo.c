#include <stdio.h>
#include <stdlib.h>

// Function to multiply two matrices using Strassen's algorithm
void strassen(int n, int A[][n], int B[][n], int C[][n]) {
    // Base case: if the matrix is 1x1, simply multiply the elements
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }

    // Divide the matrices into four submatrices
    int m = n / 2;
    int A11[m][m], A12[m][m], A21[m][m], A22[m][m];
    int B11[m][m], B12[m][m], B21[m][m], B22[m][m];
    int C11[m][m], C12[m][m], C21[m][m], C22[m][m];
    int P1[m][m], P2[m][m], P3[m][m], P4[m][m], P5[m][m], P6[m][m], P7[m][m];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
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
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            C11[i][j] = P1[i][j] + P2[i][j];
            C12[i][j] = P3[i][j] + P4[i][j];
            C21[i][j] = P5[i][j] + P6[i][j];
            C22[i][j] = P7[i][j] + C22[i][j];
        }
    }

    // Combine the submatrices into the final matrix
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            C[i][j] = C11[i][j];
            C[i][j + m] = C12[i][j];
            C[i + m][j] = C21[i][j];
            C[i + m][j + m] = C22[i][j];
        }
    }
}

// Function to print a matrix
void print_matrix(int n, int A[][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }
}

int main() {
    // Example usage
    int A[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    int B[4][4] = {{17, 18, 19, 20}, {21, 22, 23, 24}, {25, 26, 27, 28}, {29, 30, 31, 32}};
    int C[4][4];
    strassen(4, A, B, C);

    printf("A:\n");
    print_matrix(4, A);

    printf("B:\n");
    print_matrix(4, B);

    printf("C:\n");
    print_matrix(4, C);

    return 0;
}
