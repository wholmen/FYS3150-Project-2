#include <iostream>
#include <armadillo>
#include <cmath>
#include "lib.h"
#include <printf.h>
#include <time.h>
#include <string>
#include <fstream>
#include "Jacobi_function.h"

using namespace std;
using namespace arma;

mat Jacobi(mat A, int N,int *n);
void Analytical(mat A, int N);
void tqli(double *, double *, int, double **);

int main()
{

    // ----------------------------------------------------------------------
    // Implementing Jacobis method to solve an Schrodinger's equation
    // ----------------------------------------------------------------------
    int N = 1e3; int n = 0;
    //double planck = 6.62606957e-34; // m2 kg / s
    double pmax = 5;
    double h = pmax / N;
    double V;double t1;double t2;

    clock_t start1, finish1;
    start1 = clock();

    // Defining matrix A to be LHS of schrodinger equation. 
    mat A = zeros<mat>(N,N);
    A.diag(1) += -1 / (h*h); A.diag(-1) += -1 / (h*h);
    for (int i=0; i<N; i++){
        V = i*i * h*h;
        A(i,i) = 2 / (h*h) + V;
    }

    // Rotating matrix A to make it a diagonal matrix with eigenvalues along the diagonal
    A = Jacobi(A,N,&n);

    finish1 = clock();
    t1 = ( (finish1 - start1) / CLOCKS_PER_SEC );

    // ----------------------------------------------------------------
    // ------- Implementing the Solver from lib.cpp -------------------
    // ----------------------------------------------------------------

    clock_t start2, finish2;
    start2 = clock();

    // LHS of Schrodinger eq. d is the diagonal of A and e is the off-diagonal.
    double d[N]; double e[N];
    for (int i=0; i<N; i++){
        V = i*i * h*h;
        d[i] = 2 / (h*h) + V; e[i] = -1 / (h*h);
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

    finish2 = clock();
    t2 = ((finish2 - start2) /CLOCKS_PER_SEC);

    // -------------------------------------------------------------------
    // ------------ Printing results to file -----------------------------
    // -------------------------------------------------------------------

    ofstream myfile1;
    myfile1.open ("Jacobi3.txt");
    myfile1 << t1 << " " << N << " " << n << endl;
    for (int i=0;i<N-1;i++){
        myfile1 << A(i,i) << " " << A(i,i+1) << endl;
    }
    myfile1.close();

    ofstream myfile2;
    myfile2.open ("tqli3.txt");
    myfile2 << t2 << " " << N << " " << 0 << endl;
    for (int i=0;i<N;i++){
        myfile2 << d[i] << " " << e[i] << endl;
    }
    myfile2.close();

    // ------------------------------------------------------------------
    // -- Solving the Schrodinger equation for a two-particle system. ---
    // ------------------------------------------------------------------

    double wr[4]; // The oscillation fequency for the system.
    wr[0] = 0.01; wr[1] = 0.5; wr[2] = 1; wr[3] = 5;
    N = 100;
    for (int j=0;j<4;j++){ // A loop over different wr.

        // Initializing the new diagonal matrix with an updated potential depending on wr.
        double d[N]; double e[N];
        for (int i=0; i<N; i++){
            V = (i)*(i) * h*h * wr[j]*wr[j];
            d[i] = 2 / (h*h) + V;
            e[i] = -1 / (h*h);
        }

        // The eigenvector matrix needed to call function tqli
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
        ofstream myfile;

        if (j==0){myfile.open ("wr3_0,01.txt");}
        else if (j==1){myfile.open ("wr3_0,5.txt");}
        else if (j==2){myfile.open ("wr3_1.txt");}
        else {myfile.open ("wr3_5.txt");}
        myfile << 0 << " " << N << " " << 0 << endl;
        for (int k=0; k<N; k++){
            myfile << d[k] << " " << e[k] << " " << z[0][k] << " " << z[1][k] << " " << z[2][k] << endl;
        }
        myfile.close();
    }
    return 0;
}
