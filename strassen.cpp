/* strassen.cpp - an implementation of strassen's matrix multiplication algorithm */
#include <iostream>
#include <vector>
#include <limits>
using namespace std;


// S = A - B
int subtract(vector<vector<int> > &A, vector<vector<int> > &B, vector<vector<int> > &S)
{
    size_t i, j;
    if (A.size() != B.size() || 
        A[0].size() != B[0].size() || 
        A.size() != S.size() || 
        A[0].size() != S[0].size()) {
        return 0;
    }


    for (i = 0; i < A.size(); i++) {
        for (j = 0; j < B[0].size(); j++) {
            S[i][j] = A[i][j] - B[i][j];
        }
    }

    return 1;
}

// S = A + B
int add(vector<vector<int> > &A, vector<vector<int> > &B, vector<vector<int> > &S)
{
    size_t i, j;
    if (A.size() != B.size() || 
        A[0].size() != B[0].size() || 
        A.size() != S.size() || 
        A[0].size() != S[0].size()) {
        return 0;
    }

    for (i = 0; i < A.size(); i++) {
        for (j = 0; j < B[0].size(); j++) {
            S[i][j] = A[i][j] + B[i][j];
        }
    }

    return 1;
}

// C = A * B
void strassen(vector<vector<int> > &A, vector<vector<int> > &B, vector<vector<int> > &C, int n)
{
    int i, j;
    int m;

    if (n == 1) {
        C[0][0] = C[0][0] + (A[0][0] * B[0][0]);
        return;
    }

    m = n/2;

    vector<vector<int>> A11, A12, A21, A22;
    vector<vector<int>> B11, B12, B21, B22;
    vector<vector<int>> C11, C12, C21, C22;
    vector<vector<int>> S1, S2, S3, S4, S5, S6, S7, S8, S9, S10;

    // Split
    A11.resize(m, vector<int>(m));
    A12.resize(m, vector<int>(m));
    A21.resize(m, vector<int>(m));
    A22.resize(m, vector<int>(m));

    B11.resize(m, vector<int>(m));
    B12.resize(m, vector<int>(m));
    B21.resize(m, vector<int>(m));
    B22.resize(m, vector<int>(m));

    C11.resize(m, vector<int>(m));
    C12.resize(m, vector<int>(m));
    C21.resize(m, vector<int>(m));
    C22.resize(m, vector<int>(m));

    for (i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j+m];
            A21[i][j] = A[i+m][j];
            A22[i][j] = A[i+m][j+m];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j+m];
            B21[i][j] = B[i+m][j];
            B22[i][j] = B[i+m][j+m];

            C11[i][j] = C[i][j];
            C12[i][j] = C[i][j+m];
            C21[i][j] = C[i+m][j];
            C22[i][j] = C[i+m][j+m];
        }
    }

    S1.resize(m, vector<int>(m));
    S2.resize(m, vector<int>(m));
    S3.resize(m, vector<int>(m));
    S4.resize(m, vector<int>(m));
    S5.resize(m, vector<int>(m));
    S6.resize(m, vector<int>(m));
    S7.resize(m, vector<int>(m));
    S8.resize(m, vector<int>(m));
    S9.resize(m, vector<int>(m));
    S10.resize(m, vector<int>(m));

    subtract(B12, B22, S1);
    add(A11, A12, S2);
    add(A21, A22, S3);
    subtract(B21, B11, S4);
    add(A11, A22, S5);
    add(B11, B22, S6);
    subtract(A12, A22, S7);
    add(B21, B22, S8);
    subtract(A11, A21, S9);
    add(B11, B12, S10);

    vector<vector<int>> P1, P2, P3, P4, P5, P6, P7;

    P1.resize(m, vector<int>(m, 0));
    P2.resize(m, vector<int>(m, 0));
    P3.resize(m, vector<int>(m, 0));
    P4.resize(m, vector<int>(m, 0));
    P5.resize(m, vector<int>(m, 0));
    P6.resize(m, vector<int>(m, 0));
    P7.resize(m, vector<int>(m, 0));

    strassen(A11, S1, P1, m);
    strassen(S2, B22, P2, m);
    strassen(S3, B11, P3, m);
    strassen(A22, S4, P4, m);
    strassen(S5, S6, P5, m);
    strassen(S7, S8, P6, m);
    strassen(S9, S10, P7, m);

    // C11 = C11 + P5 + P4 - P2 + P6
    add(C11, P5, C11);
    add(C11, P4, C11);
    add(C11, P6, C11);
    subtract(C11, P2, C11);

    // C12 = C12 + P1 + P2
    add(C12, P1, C12);
    add(C12, P2, C12);

    // C21 = C21 + P3 + P4
    add(C21, P3, C21);
    add(C21, P4, C21);

    // C22 = C22 + P5 + P1 - P3 - P7
    add(C22, P5, C22);
    add(C22, P1, C22);
    subtract(C22, P3, C22);
    subtract(C22, P7, C22);

    // combine
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            C[i][j] = C11[i][j];
            C[i][j+m] = C12[i][j];
            C[i+m][j] = C21[i][j];
            C[i+m][j+m] = C22[i][j];
        }
    }
}


// round up to next largest power of 2
int roundUp(int n)
{
    int i;

    if (n == numeric_limits<int>::max()) return n;

    for (i = 0; i < 31; i++) {
        if ((1 << i) >= n) return (1 << i);
    }

    return -1;
}

int main()
{
    int n, p2;
    int a, b, c;
    vector<vector<int>> A, B, C;

    cout << "\nEnter Matrix Size:\n";

    cin >> n;
    if (n <= 0) {
        cout << "n must be > 0" << endl;
        return 1;
    }
    p2 = roundUp(n);

    A.resize(p2, vector<int>(p2, 0));
    B.resize(p2, vector<int>(p2, 0));
    C.resize(p2, vector<int>(p2, 0));

    cout << "\nEnter Matrix A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> a;
            A[i][j] = a;
        }
    }

    cout << "\nEnter Matrix B:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> b;
            B[i][j] = b;
        }
    }

    strassen(A, B, C, p2);

    cout << "\nA * B =\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << C[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";


    return 0;
}
