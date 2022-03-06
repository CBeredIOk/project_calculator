#include <iostream>
#include <vector>
using namespace std;
#include "Matrix_solver.h"

int main(){
    Matrix_solver Ms;

    int n, i, j, iter;
    double a, b;

    cout << "n:";
    cin >> n;

    vector<vector<double>> A(n,vector<double>(n));
    vector<double>B(n);
    vector<double>X(n);
    vector<double>Y(n);

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            A[i][j] = 0;
        }
    }
    for (j = 0; j < n; j++){
        B[j] = 0;
    }
    for (j = 0; j < n; j++){
        X[j] = 0;
    }
    for (j = 0; j < n; j++){
        Y[j] = 0;
    }

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            cout << "A(" << i << ")(" << j << ")=";
            cin >> a;
            A[i][j] = a;
        }
    }
    for (j = 0; j < n; j++){
        cout << "B(" << j << ")=";
        cin >> b;
        B[j] = b;
    }

    Ms.set_matrix(A,B,X,Y,n);

    //vector<double> ans = Ms.Gauss();
    //vector<double> ans = Ms.LU();
    //vector<double> ans = Ms.Seidel();
    vector<double> ans = Ms.Relaxation();

    //vector<vector<double>> ans_of_multiplication = Ms.Matrix_Multiplication(A,A,n);

    cout << "Answer:";
    for (j = 0; j < n; j++){
        cout << " " << ans[j] << ";";
    }

    return 0;
}