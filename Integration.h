#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

class Integration{
public:
    double f, x, y, S, a, b, p, n, k, a0, a1, a2, delta_x;
    vector<vector<double>>C, YY;

// define function

    double Function(double x){
        f = 4*x*x*x*x*x*x + 6*x*x*x*x + x;
        return f;
    }

// array X and Y

    vector<double> set_of_X(double a, double b, double p){
        delta_x = (a+b)/(2*p);
        vector<double> X(2*p+1);
        x = a;

        for(int i = 0; i <= 2*p; i++){
            X[i] = x;
            x += delta_x;
        }

        return X;
    }

    vector<double> set_of_Y(double a, double b, double p){
        delta_x = (a+b)/(2*p);
        vector<double> Y(2*p+1);
        x = a;

        for(int i = 0; i <= 2*p; i++){
            y = Function(x);
            Y[i] = y;
            x += delta_x;
        }
        return Y;
    }

// rectangle integration

    double rectangle_integration(vector<double>X,vector<double>Y){
        S = 0;
        for(int i = 1; i < X.size(); i++){
            S += (X[i]-X[i-1])*Y[i];
        }
        return S;
    }

// trapezoidal integration

    double trapezoidal_integration(vector<double>X,vector<double>Y){
        S = 0;
        for(int i = 1; i < X.size(); i++){
            S += (X[i]-X[i-1])*(Y[i]+Y[i-1])/2;
        }
        return S;
    }

// simpson integration

    double simpson_integration(vector<double>X,vector<double>Y) {
        S = 0;
        for(int i = 1; i < Y.size(); i += 2){
            S += ((X[i+1]-X[i-1])/2)*(Y[i-1]/3 + 4*Y[i]/3 + Y[i+1]/3);
        }
        return S;
    }

// lobatto integration

    double lobatto_integration(vector<double>X,vector<double>Y) {
        S = 0;

        vector<vector<double>>A_k(3,vector<double>(1));
        A_k[0][0] = 1.0/3;
        A_k[1][0] = 4.0/3;
        A_k[2][0] = 1.0/3;

        vector<vector<double>>X_k(3,vector<double>(1));
        X_k[0][0] = -1;
        X_k[1][0] = 0;
        X_k[2][0] = 1;

        vector<double> X_1(3);

        for(int i = 0; i < X.size()-1; i += 2) {
            X_1[0] = X[i];
            X_1[1] = X[i+1];
            X_1[2] = X[i+2];

            for(int j = 0; j < 3; j ++){
                S += A_k[j][0]*Function((X_1[2]-X_1[0])*X_k[j][0]/2 + (X_1[0]+X_1[2])/2);
            }
        }

        S = (X_1[2]-X_1[0])*S/2;
        return S;
    }

};