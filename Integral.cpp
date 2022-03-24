#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Integration.h"

int main(){
    Integration I;

    double a, b, err, S, S_previous, eps;
    int parts;
    vector<double> X, Y, X1, Y1;

    cout << "a:";
    cin >> a;
    cout << "b:";
    cin >> b;

    S = 0;
    parts = 2;
    err = 1;

    cout << "eps:";
    cin >> eps;

    while(err > eps){
        S_previous = S;
        X = I.set_of_X(a, b, parts);
        Y = I.set_of_Y(a, b, parts);

        //S = I.rectangle_integration(X,Y);
        //S = I.trapezoidal_integration(X,Y);
        //S = I.simpson_integration(X,Y);
        S = I.lobatto_integration(X,Y);
        
        parts = parts*3;
        err = S - S_previous;
        err = abs(err);
    }

    cout << "Answer: " << setprecision(12) << S << endl;
    return 0;
}
