#include <iostream>
#include <armadillo>
#include <cmath>
#include "lib.h"
#include <printf.h>
#include <time.h>

using namespace std;
using namespace arma;

mat Jacobi(mat A, int N);
void Analytical(mat A, int N);

int main()
{
    // Implementing Jacobis method to solve an Schrodinger's equation
    int N = 1000;
    double planck = 6.62606957e-34; // m2 kg / s
    double pmax = 1000;
    double h = pmax / N;
    double V;


    // Defining matrix A to be LHS of schrodinger equation. 

    mat A = zeros<mat>(N,N);
    A.diag(1) += -1; A.diag(-1) += -1;

    for (int i=0;i<N;i++){
        V = pow(i*h,2);
        A(i,i) = 2 + V;
    }

    A = Jacobi(A,N);
    cout << A(1,1) << endl << A(2,2) << endl << A(3,3) << endl;

    // Implementing the Solver from lib.cpp
    /*
    double d[N]; double e[N]; double *dd; double *ee;
    for (int i=0;i<N;i++){
        d[i] = 2; e[i] = -1;
    }

    dd = &d[0]; ee = &e[0];

    double **zz;
    zz = new double*[N];
    for(int i=0;i<N;i++){
        zz[i] = new double[N];
    }

    tqli(dd,ee,N,zz);
    */

    return 0;
}





mat Jacobi(mat A, int N){
    // Here the Jacobi method will be implemented. A call to this function will solve
    // an eigenvalue problem by Jacobi's method.
    double maximum = 1;
    double tolerance = 1e-4;
    int n = 0; int maxn = 1e3;


    while (maximum > tolerance && n < maxn){

        // A loop to find the highest off-diagonal element and its indices.
        int i; int j; maximum = 0.0; int I; int J;
        for (i=0; i<N; i++){
            for(j=0; j<N; j++){
                if (i != j){
                    if (abs(A(i,j)) > maximum){
                        maximum = abs(A(i,j)); I = i; J = j;
                    }
                }
            }
        } // Loop ended. Maximum element is located at (I,J)

        // Computing tau, tangens, cosine and sine for the transformation
        double tau; double t1; double t2; double t; double c; double s;
        tau = (A(J,J)-A(I,I)) / (2*A(J,I));

        // Computing tangens, and choosing the lowest of the possible roots.
        t1 = -tau + sqrt(1+pow(tau,2)); t2 = -tau - sqrt(1+pow(tau,2));
        if (abs(t1) < abs(t2)){
            t = t1;
        }
        else{
            t = t2;
        }

        // Computing cosine and sine
        c = 1/(sqrt(1+pow(t,2)));
        s = t*c;

        // Now computing the similarity matrix S and its transpose.
        mat S = eye<mat>(N,N); mat St = eye<mat>(N,N);
        S(J,J) = c; S(I,I) = c; S(J,I) = -s; S(I,J) = s;
        St(J,J) = c; St(I,I) = c; St(J,I) = s; St(I,J) = -s;

        A = St*A*S;
        n++;

    } // Ending the while-loop
    if (n == maxn){
        cout << "Jacobi did not converge fast enough. While-loop interrupted by too large n" << endl;
    }
    cout << maximum << endl;
    return A;


} // End of Jacobi function

void Analytical(mat A, int N){
    // This is a function made to compute a closed form solution of the same physical problem.

}
