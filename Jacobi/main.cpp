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
void tqli(double *, double *, int, double **);

int main()
{
    // ----------------------------------------------------------------------
    // Implementing Jacobis method to solve an Schrodinger's equation
    // ----------------------------------------------------------------------
    int N = 1000;
    //double planck = 6.62606957e-34; // m2 kg / s
    double pmax = N;
    double h = pmax / N;
    double V;

    // Defining matrix A to be LHS of schrodinger equation. 
    mat A = zeros<mat>(N,N);
    A.diag(1) += -1/h/h; A.diag(-1) += -1/h/h;
    for (int i=0;i<N;i++){
        V = (i+1)*(i+1)*h*h;
        A(i,i) = 2/h/h + V;
    }

    // Rotating matrix A to make it a diagonal matrix with eigenvalues along the diagonal
    A = Jacobi(A,N);

    // Writing out the first eigenvalues.
    cout << A(0,0) << endl << A(1,1) << endl << A(2,2) << endl << A(3,3) << endl;

    // ----------------------------------------------------------------
    // Implementing the Solver from lib.cpp
    // ----------------------------------------------------------------

    // LHS of Schrodinger eq. d is the diagonal of A and e is the off-diagonal.
    double d[N]; double e[N];
    for (int i=0;i<N;i++){
        V = i*i*h*h;
        d[i] = 2/h/h+V; e[i] = -1/h/h;
    }

    // Defining the eigenvector matrix.
    double z[N][N];
    for (int i=0;i<N;i++){
        z[i][i] = 1;
    }
    double **zz;
    zz = new double*[N];
    for(int i=0;i<N;i++){
        zz[i] = &z[i][0];
    }

    tqli(&d[0],&e[0],N,zz);
    cout << d[0] << endl << d[1] << endl << d[2];

    // -------------------------------------------------------------------
    // End of main.cpp
    // -------------------------------------------------------------------
    return 0;
}

mat Jacobi(mat A, int N){
    // Here the Jacobi method will be implemented. A call to this function will solve
    // an eigenvalue problem by Jacobi's method.
    double maximum = 1;
    double tolerance = 1e-4;
    int n = 0; int maxn = 2e5;


    while (maximum > tolerance && n < maxn){

        // A loop to find the highest off-diagonal element and its indices.
        int i; int j; maximum = 0.0; int l; int k;
        for (i=0; i<N; i++){
            for(j=i+1; j<N; j++){
                if (fabs(A(i,j)) > maximum){
                        maximum = fabs(A(i,j)); l = i; k = j;
                }
            }
        } // Loop ended. Maximum element is located at (I,J)

        // Computing tau, tangens, cosine and sine for the transformation
        double s,c;
        if (A(k,l) != 0){
            double tau; double t;
            tau = (A(l,l)-A(k,k)) / (2*A(k,l));

            // Computing tangens, and choosing the lowest of the possible roots.
            if (tau > 0){
                t = -tau + sqrt(1.0+tau*tau);
            } else {
                t = -tau - sqrt(1.0+tau*tau);
            }

            // Computing cosine and sine
            c = 1.0/sqrt(1+t*t);
            s = t*c;
        } else {
            c = 1.0; s = 0.0;
        }

        // Computing the new matrix elements.
        double a_kk, a_ll, a_ik, a_il;
        a_kk = A(k,k); a_ll = A(l,l);

        // Changing the four elements with indices l and k.
        A(k,k) = c*c*a_kk - 2.0*A(k,l)*c*s + s*s*a_ll;
        A(l,l) = s*s*a_kk + 2.0*A(k,l)*c*s + c*c*a_ll;
        A(k,l) = 0.0; A(l,k) = 0;

        // Looping over the rest of the matrix, and rotating the whole matrix.
        for (i=0; i<n; i++){
            if (i != k && i != l){
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k) = c*a_ik - s*a_il;
                A(k,i) = A(i,k);
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l);
            }
        }

        n++; // Updating n in the while-loop.

    } // Ending the while-loop
    if (n == maxn){
        cout << "Jacobi did not converge fast enough. While-loop interrupted by too large n" << endl;
    }
    return A;
} // End of Jacobi function

/*
void Analytical(mat A, int N){
    // This is a function made to compute a closed form solution of the same physical problem.

}
*/
