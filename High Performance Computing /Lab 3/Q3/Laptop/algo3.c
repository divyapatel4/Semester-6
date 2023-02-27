#include <stdio.h>
#include <stdlib.h>

#define ll long long int
#define lld long double
#define vi vector<int>
#define vd vector<double>
#define vll vector<long long int>
#define vld vector<long double>
#define vvll vector<vector<ll> >
#define vvld vector<vector<lld> >
#define vpll vector<pair<ll,ll> >
#define pll pair<ll, ll>
#define ff first
#define ss second
#define input(vec,n) for(ll i=0;i<n;i++) cin>>vec[i]
#define pb push_back
#define all(v) v.begin(), v.end()
#define endl "\n"
#define FOR(i, a, b) for (int i = (a); i < (b); ++i)
#define FOR1 (i, a) for (int i = 1; i <= (a); ++i)
#define FOR0(i, a) for (int i = 0; i < (a); ++i)
#define ROF(i, a, b) for (int i = (b); i >= (a); --i)
#define ROF0(i, a) for (int i = (b); i >= 0; --i)
#define debug(x) for(auto element:x) {cout<<element<<" ";}    cout<<endl
#define debug1(x) cout << #x << " : " << x << endl;
#define debug2(x, y) cout << #x << " : " << x << "\t" << #y << " : " << y << endl;
#define debug3(x, y, z) cout << #x << " : " << x << "\t" << #y << " : " << y << "\t" << #z << " : " << z << endl;
#define debug4(x, y, z, w) cout << #x << " : " << x << "\t" << #y << " : " << y << "\t" << #z << " : " << z << "\t" << #w << " : " << w << endl;
const ll mod = 1e9 + 7;
const ll inf = 1e18;

#define dype double

void block_matmul(int n, dype **A, dype **B, dype **C, int block_size)
{
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += block_size)
    {
        for (j = 0; j < n; j += block_size)
        {
            for (k = 0; k < n; k += block_size)
            {
                for (i1 = i; i1 < i + block_size; i1++)
                {
                    for (j1 = j; j1 < j + block_size; j1++)
                    {
                        for (k1 = k; k1 < k + block_size; k1++)
                        {
                            C[i1][j1] += A[i1][k1] * B[k1][j1];
                        }
                    }
                }
            }
        }
    }
}



int main()
{
    int n, i, j;
    scanf("%d", &n);
    dype **A = (dype **)malloc(n * sizeof(dype *));
    dype **B = (dype **)malloc(n * sizeof(dype *));
    dype **C = (dype **)malloc(n * sizeof(dype *));
    for (i = 0; i < n; i++)
    {
        A[i] = (dype *)malloc(n * sizeof(dype));
        B[i] = (dype *)malloc(n * sizeof(dype));
        C[i] = (dype *)malloc(n * sizeof(dype));
    }



    // matA = np.array([ [ 1, 2, 3, 4 ], [ 5, 6, 7, 8 ], [ 9, 10, 11, 12 ], [ 13, 14, 15, 16 ] ])
    // matB = np.array([ [ 17, 18, 19, 20 ], [ 21, 22, 23, 24 ], [ 25, 26, 27, 28 ], [ 29, 30, 31, 32 ] ])

    A[0][0] = 1; A[0][1] = 2; A[0][2] = 3; A[0][3] = 4; A[1][0] = 5; A[1][1] = 6; A[1][2] = 7; A[1][3] = 8; A[2][0] = 9; A[2][1] = 10; A[2][2] = 11; A[2][3] = 12; A[3][0] = 13; A[3][1] = 14; A[3][2] = 15; A[3][3] = 16;
    B[0][0] = 17; B[0][1] = 18; B[0][2] = 19; B[0][3] = 20; B[1][0] = 21; B[1][1] = 22; B[1][2] = 23; B[1][3] = 24; B[2][0] = 25; B[2][1] = 26; B[2][2] = 27; B[2][3] = 28; B[3][0] = 29; B[3][1] = 30; B[3][2] = 31; B[3][3] = 32;

    block_matmul(n, A, B, C, 2);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d ", C[i][j]);
        }
        printf("\n");
    }

}
