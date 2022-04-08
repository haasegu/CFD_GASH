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



void JacobiSolve(CRS_Matrix1 const &SK1, vector<double> const &f, vector<double> &u)
{
    const double omega   = 1.0;
    const int    maxiter = 100;
    const double tol  = 1e-5,                // tolerance
                 tol2 = tol * tol;           // tolerance^2

    int nrows=SK1.Nrows();                    // number of rows == number of columns
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

    SK1.GetDiag(dd);                          //  dd := diag(K)
    ////DebugVector(dd);{int ijk; cin >> ijk;}

    //  Initial sweep
    SK1.Defect(r, f, u);                      //  r := f - K*u

    vddiv(w, r, dd);                         //  w := D^{-1}*r
    double sigma0 = dscapr(w, r);            // s0 := <w,r>

    // Iteration sweeps
    int iter  = 0;
    double sigma = sigma0;
    while ( sigma > tol2 * sigma0 && maxiter > iter)
    {
        ++iter;
        vdaxpy(u, u, omega, w );             //  u := u + om*w
        SK1.Defect(r, f, u);                  //  r := f - K*u
        vddiv(w, r, dd);                     //  w := D^{-1}*r
        sigma = dscapr(w, r);                // s0 := <w,r>
      	cout << "Iteration " << iter << " : " << sqrt(sigma/sigma0) << endl;
    }
    cout << "aver. Jacobi rate :  " << exp(log(sqrt(sigma / sigma0)) / iter) << "  (" << iter << " iter)" << endl;
    cout << "final error: " << sqrt(sigma / sigma0) << " (rel)   " << sqrt(sigma) << " (abs)\n";

    return;
}








void JacobiSmoother(Matrix const &SK1, std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> &r, int nsmooth, double const omega, bool zero)
{
    //// ToDO: ensure compatible dimensions
    SK1.JacobiSmoother(f, u, r, nsmooth, omega, zero);
    return;

    int const nnodes = static_cast<int>(u.size());
    if (zero) {            // assumes initial solution is zero
        DiagPrecond(SK1, f, u, omega);
        --nsmooth;                           // first smoothing sweep done
    }
    //cout << zero << endl;

    auto const &D = SK1.GetDiag();            // accumulated diagonal of matrix @p SK.
    //auto const D = SK.GetDiag();            // accumulated diagonal of matrix @p SK.
    for (int ns = 1; ns <= nsmooth; ++ns) {
        SK1.Defect(r, f, u);                  //  r := f - K*u
#pragma omp parallel for    
        for (int k = 0; k < nnodes; ++k) {
            // u := u + om*D^{-1}*r
            u[k] = u[k] + omega * r[k] / D[k]; // MPI: distributed to accumulated vector needed
        }
    }

    return;
}

void DiagPrecond(Matrix const &SK1, std::vector<double> const &r, std::vector<double> &w,
                 double const omega)
{
    // ToDO: ensure compatible dimensions
    auto const &D = SK1.GetDiag();        // accumulated diagonal of matrix @p SK.
    int const nnodes = static_cast<int>(w.size());
#pragma omp parallel for    
    for (int k = 0; k < nnodes; ++k) {
        w[k] = omega * r[k] / D[k];      // MPI: distributed to accumulated vector needed
    }

    return;
}



