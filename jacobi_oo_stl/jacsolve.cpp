#include "vdop.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "userset.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <colamd.h>
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



/*

Multigrid::Multigrid(Mesh const &cmesh, int const nlevel)
    : _meshes(cmesh, nlevel),
      _SK1(), _u(_meshes.size()), _f(_meshes.size()), _d(_meshes.size()), _w(_meshes.size()),
      _Pc2f()
{
    cout << "\n........................  in Multigrid::Multigrid  ..................\n";
    // Allocate Memory for matrices/vectors on all levels
    for (size_t lev = 0; lev < Nlevels(); ++lev) {
        _SK1.push_back( FEM_Matrix(_meshes[lev]) );  // CRS matrix
        const auto nn = _SK1[lev].Nrows();
        _u[lev].resize(nn);
        _f[lev].resize(nn);
        _d[lev].resize(nn);
        _w[lev].resize(nn);
        auto vv = _meshes[lev].GetFathersOfVertices();
        cout << vv.size() << endl;
    }
    // Intergrid transfer operators
    //cout << "\n........................  in Multigrid::Multigrid  Prolongation ..................\n";
    //_Pc2f.push_back( BisectInterpolation(vector<int>(0)) ); // no prolongation to coarsest grid
    _Pc2f.push_back( BisectIntDirichlet() ); // no prolongation to coarsest grid
    for (size_t lev = 1; lev < Nlevels(); ++lev) {
                    //cout << lev << endl;
                    //cout << _meshes[lev].GetFathersOfVertices () << endl;
        _Pc2f.push_back( BisectIntDirichlet( _meshes[lev].GetFathersOfVertices (), _meshes[lev-1].Index_DirichletNodes ()  )  );
                    //cout << _Pc2f.back().Nrows() << "  " << _Pc2f.back().Ncols() << endl;
    }
    cout << "\n..........................................\n";
}

Multigrid::~Multigrid()
{}

void Multigrid::DefineOperators()
{
    for (size_t lev = 0; lev < Nlevels(); ++lev) {
        DefineOperator(lev);
    }
    return;
}

// GH: Hack
void Multigrid::DefineOperator(size_t lev)
{
    _SK1[lev].CalculateLaplace(_f[lev]);  // fNice()  in userset.h

    if (lev == Nlevels() - 1) {                // fine mesh
        _meshes[lev].SetValues(_u[lev], [](double x, double y, double z) -> double
        { return x *z * std::sin(2.5 * M_PI * y); }
                              );
    }
    else 
    {
        _meshes[lev].SetValues(_u[lev], f_zero);
    }

    _SK1[lev].ApplyDirichletBC(_u[lev], _f[lev]);

    return;
}


void Multigrid::DefineOperators_heat_equation(const vector<double> &u_old, double const dt, double const t, const double c)
{
    _SK1.clear();

    for (size_t lev = 0; lev < Nlevels(); ++lev) {
        _SK1.push_back( FEM_Matrix(_meshes[lev]) );  // CRS matrix
        _meshes[lev].SetValues(_u[lev], f_zero);
        _meshes[lev].SetValues(_f[lev], f_zero);
        _meshes[lev].SetValues(_d[lev], f_zero);
        _meshes[lev].SetValues(_w[lev], f_zero);
    }

    for (size_t lev = 0; lev < Nlevels(); ++lev) {
        DefineOperator_heat_equation(lev, u_old, dt, t, c);
    }
    return;
}

void Multigrid::DefineOperator_heat_equation(size_t lev, const vector<double> &u_old, double const dt, double const t, const double c)
{
    

    _SK1[lev].CalculateLaplace_heat_equation(_f[lev], u_old, dt, t, c);  // fNice()  in userset.h

    if (lev == Nlevels() - 1) {                // fine mesh
        _meshes[lev].SetValues(_u[lev], [](double x, double y, double z) -> double
        { return 0.0*x *z * std::sin(2.5 * M_PI * y); }
                              );
    }
    else {
        _meshes[lev].SetValues(_u[lev], f_zero);
    }

    auto const &xc   = _meshes[lev].GetCoords();
    //_SK[lev].ApplyDirichletBC(_u[lev], _f[lev]);
    _SK1[lev].ApplyDirichletBC_Box(_u[lev],_f[lev],1.0,1.0,1.0,1.0,1.0,1.0);


    return;
}

void Multigrid::JacobiSolve(size_t lev)
{
    assert(lev < Nlevels());
    ::JacobiSolve(_SK1[lev], _f[lev], _u[lev]);
}

void Multigrid::MG_Step(size_t lev, int const pre_smooth, bool const bzero, int nu)
{
    assert(lev < Nlevels());
    int const post_smooth = pre_smooth;

    if (lev == 0) { // coarse level
        JacobiSmoother(_SK1[lev], _f[lev], _u[lev], _d[lev],  100, 1.0, false);
    }
    else {
        JacobiSmoother(_SK1[lev], _f[lev], _u[lev], _d[lev],  pre_smooth, 0.85, bzero);

        if (nu > 0) {

            _SK1[lev].Defect(_d[lev], _f[lev], _u[lev]);   //   d := f - K*u
            _Pc2f[lev].MultT(_d[lev], _f[lev - 1]);       // f_H := R*d
            //DefectRestrict(_SK[lev], _Pc2f[lev], _f[lev - 1], _f[lev], _u[lev]); // f_H := R*(f - K*u)

                    //_meshes[lev-1].Visualize(_f[lev - 1]);        // GH: Visualize: f_H should be 0 on Dirichlet B.C.

            MG_Step(lev - 1, pre_smooth, true, nu);       // solve  K_H * u_H =f_H  with u_H:=0
            for (int k = 1; k < nu; ++k) {
                // W-cycle
                MG_Step(lev - 1, pre_smooth, false, nu);  // solve  K_H * u_H =f_H
            }

            _Pc2f[lev].Mult(_w[lev], _u[lev - 1]);        // w := P*u_H

            vdaxpy(_u[lev], _u[lev], 1.0, _w[lev] );      // u := u + tau*w
        }

        JacobiSmoother(_SK1[lev], _f[lev], _u[lev], _d[lev],  post_smooth, 0.85, false);

    }

    return;
}

void Multigrid::MG_Solve(int pre_smooth, double eps, const bool paraview, const int iteration, int nu)
{
    size_t lev=Nlevels()-1;                // fine level

    // start with zero guess
    DiagPrecond(_SK1[lev], _f[lev], _w[lev], 1.0);  // w   := D^{-1]*f
    //double s0 = L2_scapr(_f[lev],_w[lev]);         // s_0 := <f,w>
    double s0 = dscapr(_f[lev],_w[lev]);         // s_0 := <f,w>
    double si;

    bool bzero = true;                       // start with zero guess
    int  iter  = 0;
    do
    {
        MG_Step(lev, pre_smooth, bzero, nu);
        bzero=false;
        _SK1[lev].Defect(_d[lev], _f[lev], _u[lev]);    //   d := f - K*u
        DiagPrecond(_SK1[lev], _d[lev], _w[lev], 1.0);  // w   := D^{-1]*d
        //si = L2_scapr(_d[lev],_w[lev]);                // s_i := <d,w>
        si = dscapr(_d[lev],_w[lev]);                // s_i := <d,w>
        ++iter;
    } while (si>s0*eps*eps);



    // visualize defect
    //size_t level=Nlevels()-1;
    //const auto &ml=GetMesh(level);
    //auto s = std::to_string(iteration);

    //string fname;

    // if(iteration!=-1)
    // {
    //     if(paraview==1)
    //     {
    //         fname="defect"+s+".vtk";
    //         ml.Write_ascii_paraview(fname,_d[level]);
    //     }
    //     else
    //     {
    //         fname="defect"+s+".vtk";
    //         ml.Write_ascii_matlab(fname,_d[level]);
    //     }
    // }

    cout << "\nrel. error: " << sqrt(si/s0) << "  ( " << iter << " iter.)" << endl;
    return;
}


*/



