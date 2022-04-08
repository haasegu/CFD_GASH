#include "vdop.h"
#include "geom.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "userset.h"

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
    //const int    maxiter = 1000;
    const int    maxiter = 240; // GH
    const double tol  = 1e-6,                // tolerance
                 tol2 = tol * tol;           // tolerance^2

    int nrows = SK.Nrows();                  // number of rows == number of columns
    assert( nrows == static_cast<int>(f.size()) && f.size() == u.size() );

    cout << endl << " Start Jacobi solver for " << nrows << " d.o.f.s"  << endl;
    //  Choose initial guess
    for (int k = 0; k < nrows; ++k) {
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
    const double sigma0 = dscapr(w, r);      // s0 := <w,r>

    // Iteration sweeps
    int iter  = 0;
    double sigma = sigma0;
    double sig_old=sigma;
    while ( sigma > tol2 * sigma0 && maxiter > iter)  // relative error
        //while ( sigma > tol2 && maxiter > iter)         // absolute error
    {
        sig_old=sigma;
        ++iter;
        vdaxpy(u, u, omega, w );             //  u := u + om*w
        SK.Defect(r, f, u);                  //  r := f - K*u
        vddiv(w, r, dd);                     //  w := D^{-1}*r
        sigma = dscapr(w, r);                // s0 := <w,r>
      	//cout << "Iteration " << iter << " : " << sqrt(sigma/sigma0) << endl;
        if (sigma>sig_old) cout << "Divergent at iter " << iter << endl;    // GH
    }
    cout << "aver. Jacobi rate :  " << exp(log(sqrt(sigma / sigma0)) / iter) << "  (" << iter << " iter)" << endl;
    cout << "final error: " << sqrt(sigma / sigma0) << " (rel)   " << sqrt(sigma) << " (abs)\n";

    return;
}



void JacobiSmoother(Matrix const &SK, std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> &r, int nsmooth, double const omega, bool zero)
{
    //// ToDO: ensure compatible dimensions
    SK.JacobiSmoother(f, u, r, nsmooth, omega, zero);
    return;

    int const nnodes = static_cast<int>(u.size());
    if (zero) {            // assumes initial solution is zero
        DiagPrecond(SK, f, u, omega);
        --nsmooth;                           // first smoothing sweep done
    }
    //cout << zero << endl;

    auto const &D = SK.GetDiag();            // accumulated diagonal of matrix @p SK.
    //auto const D = SK.GetDiag();            // accumulated diagonal of matrix @p SK.
    for (int ns = 1; ns <= nsmooth; ++ns) {
        SK.Defect(r, f, u);                  //  r := f - K*u
#pragma omp parallel for    
        for (int k = 0; k < nnodes; ++k) {
            // u := u + om*D^{-1}*r
            u[k] = u[k] + omega * r[k] / D[k]; // MPI: distributed to accumulated vector needed
        }
    }

    return;
}

void DiagPrecond(Matrix const &SK, std::vector<double> const &r, std::vector<double> &w,
                 double const omega)
{
    // ToDO: ensure compatible dimensions
    auto const &D = SK.GetDiag();        // accumulated diagonal of matrix @p SK.
    int const nnodes = static_cast<int>(w.size());
#pragma omp parallel for    
    for (int k = 0; k < nnodes; ++k) {
        w[k] = omega * r[k] / D[k];      // MPI: distributed to accumulated vector needed
    }

    return;
}


Multigrid::Multigrid(Mesh const &cmesh, int const nlevel)
    : _meshes(cmesh, nlevel),
      _SK(), _u(_meshes.size()), _f(_meshes.size()), _d(_meshes.size()), _w(_meshes.size()),
      _Pc2f()
{
    cout << "\n........................  in Multigrid::Multigrid  ..................\n";
    // Allocate Memory for matrices/vectors on all levels
    for (size_t lev = 0; lev < Nlevels(); ++lev) {
        _SK.push_back( FEM_Matrix(_meshes[lev]) );  // CRS matrix
        const auto nn = _SK[lev].Nrows();
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
        //_Pc2f.push_back( BisectInterpolation( _meshes[lev].GetFathersOfVertices () )  );
        _Pc2f.push_back( BisectIntDirichlet( _meshes[lev].GetFathersOfVertices (), _meshes[lev-1].Index_DirichletNodes ()  )  );
                    //cout << _Pc2f.back().Nrows() << "  " << _Pc2f.back().Ncols() << endl;
        //checkInterpolation(lev);
        //checkRestriction(lev);
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
#include "omp.h"
void Multigrid::DefineOperator(size_t lev)
{
    double tstart = omp_get_wtime(); 
    _SK[lev].CalculateLaplace(_f[lev]);  // fNice()  in userset.h
    double t1 = omp_get_wtime() - tstart;             // OpenMP
    cout << "CalculateLaplace: timing in sec. : " << t1 << "   in level " << lev << endl;

    if (lev == Nlevels() - 1) {                // fine mesh
        _meshes[lev].SetValues(_u[lev], [](double x, double y) -> double
        { return x *x * std::sin(2.5 * M_PI * y); }
                              );
    }
    else {
        _meshes[lev].SetValues(_u[lev], f_zero);
    }

//// GH:
    ////_SK[lev].CheckRowSum();
    ////if (!_SK[lev].CheckMproperty()) _SK[lev].ForceMproperty();
    
    //assert(_SK[lev].CheckSymmetry());
    //assert(_SK[lev].CheckRowSum());
    ////assert(_SK[lev].CheckMproperty());
    //_SK[lev].CheckMproperty();
//// HG
    
    _SK[lev].ApplyDirichletBC(_u[lev], _f[lev]);

    return;
}

void Multigrid::JacobiSolve(size_t lev)
{
    assert(lev < Nlevels());
    ::JacobiSolve(_SK[lev], _f[lev], _u[lev]);
}

void Multigrid::MG_Step(size_t lev, int const pre_smooth, bool const bzero, int nu)
{
    assert(lev < Nlevels());
    int const post_smooth = pre_smooth;

    if (lev == 0) { // coarse level
        // GH: a factorization (once in setup) with repeated forward-backward substitution would be better
        int n_jacobi_iterations = GetMesh(lev).Nnodes()/10; // ensure accuracy for coarse grid silver
        JacobiSmoother(_SK[lev], _f[lev], _u[lev], _d[lev],  n_jacobi_iterations, 1.0, false);
        //JacobiSmoother(_SK[lev], _f[lev], _u[lev], _d[lev],  1000, 1.0, false);
    }
    else {
        JacobiSmoother(_SK[lev], _f[lev], _u[lev], _d[lev],  pre_smooth, 1.0, bzero || lev < Nlevels()-1);
        //JacobiSmoother(_SK[lev], _f[lev], _u[lev], _d[lev],  pre_smooth, 0.85, bzero);

        if (nu > 0) {

            _SK[lev].Defect(_d[lev], _f[lev], _u[lev]);   //   d := f - K*u
            _Pc2f[lev].MultT(_d[lev], _f[lev - 1]);       // f_H := R*d
            //DefectRestrict(_SK[lev], _Pc2f[lev], _f[lev - 1], _f[lev], _u[lev]); // f_H := R*(f - K*u)
            
            //cout << "f fine  " << endl; GetMesh(lev).Visualize(_f[lev]);
            //cout << "u fine  " << endl; GetMesh(lev).Visualize(_u[lev]);
            //cout << "d fine  " << endl; GetMesh(lev).Visualize(_d[lev]);
            //cout << "f coarse" << endl; GetMesh(lev-1).Visualize(_f[lev - 1]);

                    //_meshes[lev-1].Visualize(_f[lev - 1]);        // GH: Visualize: f_H should be 0 on Dirichlet B.C.

            MG_Step(lev - 1, pre_smooth, true, nu);       // solve  K_H * u_H =f_H  with u_H:=0
            for (int k = 1; k < nu; ++k) {
                // W-cycle
                MG_Step(lev - 1, pre_smooth, false, nu);  // solve  K_H * u_H =f_H
            }

            _Pc2f[lev].Mult(_w[lev], _u[lev - 1]);        // w := P*u_H

            vdaxpy(_u[lev], _u[lev], 1.0, _w[lev] );      // u := u + tau*w
        }

        JacobiSmoother(_SK[lev], _f[lev], _u[lev], _d[lev],  post_smooth, 0.85, false);
        
        //// Check Interpolation
        //vector<double> const vc(_u[lev - 1].size(),1.0);
        //vector<double>       vf(_u[lev].size(),-1.0);
        //_Pc2f[lev].Mult(vf, vc);  
        //GetMesh(lev).Visualize(vf);

        //// Check Restriction
        //vector<double>       vc(_u[lev - 1].size(),-1.0);
        //vector<double> const vf(_u[lev].size(),1.0);
        //_Pc2f[lev].MultT(vf, vc);  
        //GetMesh(lev-1).Visualize(vc);


    }

    return;
}

void Multigrid::MG_Solve(int pre_smooth, double eps, int nu)
{
    size_t lev=Nlevels()-1;                // fine level

    // start with zero guess
    DiagPrecond(_SK[lev], _f[lev], _w[lev], 1.0);  // w   := D^{-1]*f
    _u[lev] = _w[lev];
    
    //double s0 = L2_scapr(_f[lev],_w[lev]);         // s_0 := <f,w>
    double s0 = dscapr(_f[lev],_w[lev]);         // s_0 := <f,w>
    double si;

    bool bzero = true;                       // start with zero guess
    int  iter  = 0;
    do
    {
        MG_Step(lev, pre_smooth, bzero, nu);
        bzero=false;
        _SK[lev].Defect(_d[lev], _f[lev], _u[lev]);    //   d := f - K*u
        DiagPrecond(_SK[lev], _d[lev], _w[lev], 1.0);  // w   := D^{-1]*d
        //si = L2_scapr(_d[lev],_w[lev]);                // s_i := <d,w>
        si = dscapr(_d[lev],_w[lev]);                // s_i := <d,w>
        ++iter;
    } while (si>s0*eps*eps);


    cout << "\nrel. error: " << sqrt(si/s0) << "  ( " << iter << " iter.)" << endl;
    return;
}

bool Multigrid::checkInterpolation(size_t const lev)
{
    assert(1<=lev && lev<Nlevels());
    _meshes[lev-1].SetValues(_w[lev-1], [](double x, double y) -> double
                           { return x+y; }  );
    _meshes[lev].SetValues(_w[lev], [](double x, double y) -> double
                           { return -123.0; }  );
    //static_cast<BisectInterpolation>(_Pc2f[lev]).Mult(_d[lev], _w[lev - 1]);        // d := P*w_H
    _Pc2f[lev].Mult(_d[lev], _w[lev - 1]);        // d := P*w_H
    
    cout << "înterpolated  " << endl; GetMesh(lev).Visualize(_d[lev]);
    
    return true;
}

bool Multigrid::checkRestriction(size_t const lev)
{
    assert(1<=lev && lev<Nlevels());
    _meshes[lev].SetValues(_d[lev], [](double x, double y) -> double
                           { return x+y; }  );
    _meshes[lev-1].SetValues(_w[lev-1], [](double x, double y) -> double
                           { return -123.0; }  );
    //static_cast<BisectInterpolation>(_Pc2f[lev]).MultT(_d[lev], _w[lev - 1]);        // w_H := R*d
    _Pc2f[lev].MultT(_d[lev], _w[lev - 1]);        // w_H := R*d
    
    cout << "restricted  " << endl; GetMesh(lev-1).Visualize(_w[lev-1]);
    
    return true;
}


