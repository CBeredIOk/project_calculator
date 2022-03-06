#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

class Matrix_solver{
public:
    vector<vector<double>>A;
    vector<double>B;
    vector<double>X;
    vector<double>Y;
    double a, d, s, b, w, err, eps;
    int n, i, j, k;

// Функция выводящая матрицу в ком. строку

    void show(vector <vector <double>> A, int n){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                cout << "\t" << A[i][j] << "\t";
            }
            cout << endl;
        }
    }

// Установка матрицы

    void set_matrix(vector<vector<double>>A, vector<double>B,vector<double>X, vector<double>Y, int n){
        this->A=A;
        this->B=B;
        this->X=X;
        this->Y=Y;
        this->n=n;
        cout << "Matrix is set" << endl;
        show(A,n);
    }

// Решение СЛАУ Гауссом

    vector<double> Gauss(vector<vector<double>>A, vector<double>B){
        for (k = 0; k < n; k++){ // прямой ход
            for (j = k+1; j < n; j++){
                d = A[j][k] / A[k][k]; // формула (1)
                for (i = k; i < n; i++){
                    A[j][i] = A[j][i] - d * A[k][i]; // формула (2)
                }
                B[j] = B[j] - d * B[k]; // формула (3)
            }
        }

        for (k = n-1; k >= 0; k--){ // обратный ход
            d = 0;
            for (j = k; j < n; j++){
                s = A[k][j] * X[j]; // формула (4)
                d = d + s; // формула (4)
            }
            X[k] = (B[k] - d) / A[k][k]; // формула (4)
        }
        return X;
    };

// LU решение СЛАУ

    void LU_1(vector <vector <double>> A, vector <vector <double>> &L, vector <vector <double>> &U, int n){
        //находим первый столбец L[][] и первую строку U[][]
        double sum;
        for (int i = 0; i < n; i++){
            L[i][0] = A[i][0];
            U[0][i] = A[0][i]/L[0][0];
        }

        //дальше вычисляем L[][], U[][] по формуле

        for (int i = 1; i < n; i++){
            for (int j = 1; j < n; j++){
                if (i >= j){ //нижний треугольник
                    sum = 0;
                    for (int k = 0; k < j; k++)
                        sum += L[i][k] * U[k][j];

                    L[i][j] = A[i][j] - sum;
                }
                else{ // верхний
                    sum = 0;
                    for (int k = 0; k < i; k++)
                        sum += L[i][k] * U[k][j];

                    U[i][j] = (A[i][j] - sum) / L[i][i];
                }
            }
        }
    }

    vector<double>LU(){
        vector <vector<double>> L(n,vector<double>(n)), U(n,vector<double>(n));
        for (int i  = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                L[i][j] = 0;
                U[i][j] = 0;

                if (i == j)
                    U[i][j] = 1;
            }
        }

        LU_1(A,L,U, n);

        Y = Gauss(L,B);

        for (j = 0; j < n; j++){
            X[j] = 0;
        }

        for (k = n-1; k >= 0; k--){ // обратный ход
            d = 0;
            for (j = k; j < n; j++){
                s = U[k][j] * X[j]; // формула (4)
                d = d + s; // формула (4)
            }
            X[k] = (Y[k] - d) / U[k][k]; // формула (4)
        }

        return X;
    };

// Зейдель

    vector<double>Seidel(){
        cout << "eps:";
        cin >> eps;

        vector<double>X0(n);
        for (j = 0; j < n; j++){
            X0[j] = 0;
        }
        vector<double>X0_2(n);
        for (j = 0; j < n; j++){
            X0_2[j] = 0;
        }

        err = 1;

        while(true){
            for (i = 0; i < n; i++){
                double c = B[i]/A[i][i];
                for (j = 0; j < n; j++){
                    b = -A[i][j] / A[i][i];
                    X[i] += b * X0[j];
                }
                X[i] += c;
                X0[i] = X[i];
            }

            err = secondNorm(X,X0_2,n);

            X0_2 = X;

            err = sqrt(err);

            if(err < eps){
                break;
            }
        }
        return X;
    };

// Релаксации

    vector<double>Relaxation(){
        cout << "eps:";
        cin >> eps;
        cout << "w:";
        cin >> w;

        vector<double>X0(n);
        for (j = 0; j < n; j++){
            X0[j] = 0;
        }
        vector<double>X0_2(n);
        for (j = 0; j < n; j++){
            X0_2[j] = 0;
        }

        vector<vector<double>>HELP(n,vector<double>(n));

        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                HELP[i][j] = -w;
            }
            HELP[i][i] += 1;
        }

        err = 1;

        while(true){
            for (i = 0; i < n; i++){
                a = 0;
                for (j = 0; j < n; j++){
                    a = a + HELP[i][j]*A[i][j]/A[i][i]*X[j];
                }
                X[i] = a + w*B[i]/A[i][i];
            }

            err = secondNorm(X,X0_2,n);

            X0_2 = X;

            err = sqrt(err);

            if(err < eps){
                break;
            }
        }
        return X;
    };

// Умножение матриц

    vector<vector<double>> Matrix_Multiplication(vector<vector<double>>A, vector<vector<double>>B, int n){
        vector<vector<double>>C(n,vector<double>(n));
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                C[i][j] = 0;
            }
        }

        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                C[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    };

// Вторая норма вектора

    double secondNorm(vector<double> X,vector<double> X0, int n){
        double abs;
        for(i = 0; i < n; i++){
            abs += (X[i]-X0[i])*(X[i]-X0[i]);
        }
        abs = sqrt(abs);
        return abs;
    };

};

