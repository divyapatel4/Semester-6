#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long lld;

/* Strassen's Algorithm for matrix multiplication
Complexity: O(n^2.808) */

inline lld** MatrixMultiply(lld** a, lld** b, int n,
                            int l, int m)
{
    lld** c = (lld**)malloc(n * sizeof(lld*));
    for (int i = 0; i < n; i++)
        c[i] = (lld*)malloc(m * sizeof(lld));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            c[i][j] = 0;
            for (int k = 0; k < l; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

inline lld** Strassen(lld** a, lld** b, int n,
                      int l, int m)
{
    if (n == 1 || l == 1 || m == 1)
        return MatrixMultiply(a, b, n, l, m);

    lld** c = (lld**)malloc(n * sizeof(lld*));
    for (int i = 0; i < n; i++)
        c[i] = (lld*)malloc(m * sizeof(lld));

    int adjN = (n >> 1) + (n & 1);
    int adjL = (l >> 1) + (l & 1);
    int adjM = (m >> 1) + (m & 1);

    lld**** As = (lld****)malloc(2 * sizeof(lld***));
    for (int x = 0; x < 2; x++) {
        As[x] = (lld***)malloc(2 * sizeof(lld**));
        for (int y = 0; y < 2; y++) {
            As[x][y] = (lld**)malloc(adjN * sizeof(lld*));
            for (int i = 0; i < adjN; i++) {
                As[x][y][i] = (lld*)malloc(adjL * sizeof(lld));
                for (int j = 0; j < adjL; j++) {
                    int I = i + (x & 1) * adjN;
                    int J = j + (y & 1) * adjL;
                    As[x][y][i][j] = (I < n && J < l) ? a[I][J] : 0;
                }
            }
        }
    }

    lld**** Bs = (lld****)malloc(2 * sizeof(lld***));
    for (int x = 0; x < 2; x++) {
        Bs[x] = (lld***)malloc(2 * sizeof(lld**));
        for (int y = 0; y < 2; y++) {
            Bs[x][y] = (lld**)malloc(adjL * sizeof(lld*));
            lld ***p1 = Strassen(As[0][0], s[0], adjN, adjL, adjM);
            lld ***p2 = Strassen(As[1][1], s[1], adjN, adjL, adjM);
            lld ***p3 = Strassen(As[0][1], Bs[1][0], adjN, adjL, adjM);
            lld ***p4 = Strassen(As[1][0], Bs[0][0], adjN, adjL, adjM);
            lld ***p5 = Strassen(As[0][0], s[2], adjN, adjL, adjM);
            lld ***p6 = Strassen(As[1][1], s[3], adjN, adjL, adjM);
            lld ***p7 = Strassen(As[0][1], Bs[1][1], adjN, adjL, adjM);


            for (int i = 0; i < adjN; i++) {
                for (int j = 0; j < adjM; j++) {
                    c[i][j] = p4[i][j] + p5[i][j] + p6[i][j] - p2[i][j];
                }
            }
            for (int i = 0; i < adjN; i++) {
                for (int j = adjM; j < m; j++) {
                    c[i][j] = p1[i][j - adjM] + p2[i][j - adjM];
                }
            }
            for (int i = adjN; i < n; i++) {
                for (int j = 0; j < adjM; j++) {
                    c[i][j] = p3[i - adjN][j] + p4[i - adjN][j];
                }
            }
            for (int i = adjN; i < n; i++) {
                for (int j = adjM; j < m; j++) {
                    c[i][j] = p5[i - adjN][j - adjM] + p1[i - adjN][j - adjM] -
                                    p3[i - adjN][j - adjM] - p7[i - adjN][j - adjM];
                }
            }

return c;
}
    }
}

int main()
{
	int n, m;
	cin >> n >> m;

	lld** a = new lld*[n];
	for (int i = 0; i < n; i++)
		a[i] = new lld[m];

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cin >> a[i][j];
		}
	}

	lld** b = new lld*[n];
	for (int i = 0; i < n; i++)
		b[i] = new lld[m];

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cin >> b[i][j];
		}
