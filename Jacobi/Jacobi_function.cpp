# include "Jacobi_function.h"

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
