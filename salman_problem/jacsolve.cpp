#include "vdop.h"
#include "getmatrix.h"
#include "jacsolve.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

// #####################################################################
// 	const int neigh[], const int color,
// 	const MPI::Intracomm& icomm,
void JacobiSolve(CRS_Matrix const &SK, vector<double> const &f, vector<double> &u)
{
    const double omega   = 1.0;
    const int    maxiter = 100;
    const double tol  = 1e-5,                // tolerance
                 tol2 = tol * tol;           // tolerance^2

    int nrows=SK.Nrows();                    // number of rows == number of columns
    assert( nrows==static_cast<int>(f.size()) && f.size()==u.size() );

    cout << endl << " Start Jacobi solver for " << nrows << " d.o.f.s"  << endl;
    //  Choose initial guess
    for (int k = 0; k < nrows; ++k)
    {
        u[k] = 0.0;                          //  u := 0
    }

    vector<double> dd(nrows);                // matrix diagonal
    vector<double>  r(nrows);                // residual
    vector<double>  w(nrows);                // correction

    SK.GetDiag(dd);                          //  dd := diag(K)
    ////DebugVector(dd);{int ijk; cin >> ijk;}

    //  Initial sweep
    SK.Defect(r, f, u);                      //  r := f - K*u

    vddiv(w, r, dd);                         //  w := D^{-1}*r
    double sigma0 = dscapr(w, r);            // s0 := <w,r>

    // Iteration sweeps
    int iter  = 0;
    double sigma = sigma0;
    while ( sigma > tol2 * sigma0 && maxiter > iter)
    {
        ++iter;
        vdaxpy(u, u, omega, w );             //  u := u + om*w
        SK.Defect(r, f, u);                  //  r := f - K*u
        vddiv(w, r, dd);                     //  w := D^{-1}*r
        sigma = dscapr(w, r);                // s0 := <w,r>
      	cout << "Iteration " << iter << " : " << sqrt(sigma/sigma0) << endl;
    }
    cout << "aver. Jacobi rate :  " << exp(log(sqrt(sigma / sigma0)) / iter) << "  (" << iter << " iter)" << endl;
    cout << "final error: " << sqrt(sigma / sigma0) << " (rel)   " << sqrt(sigma) << " (abs)\n";

    return;
}

