//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "elements.h"                   // GH: finite elements collection
#include "geom.h"
//#include "geom3.h"
#include "getmatrix.h"
//#include "getmatrix_2.h"
#include "jacsolve.h"
#include "userset.h"             // rhs_lap3(), sol_lap3()
#include "vdop.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <omp.h>
#include <chrono>           // timing
using namespace std;
using namespace std::chrono;  // timing

/**
 * Checks correctness for 3D Laplace (in a box domain) 
 * using tetrahedrons with linear test function, scalar values.
 * 
 * time dependent iteration ?
 * 
 * @see P1_3d
 */
bool check_P1_3d();

/**
 * Checks correctness for 3D xxx (in a box domain) 
 * using tetrahedrons with linear test function for scalar values
 * and quadratic test functions for the vector values
 * 
 * @see P1_2vec_3d
 */
bool check_P1_2vec_3d();


int main(int, char **)
{
    //check_P1_3d();
    
    check_P1_2vec_3d();

    return 0;
}


bool check_P1_3d()
{     
    //Mesh_3d_P1_matlab const mesh("../salman_problem/NS_3D.txt");  // old
    //Mesh const mesh("../salman_problem/NS_3D.txt");               // new
    Mesh const mesh("../salman_problem/GH_NS_3D.txt");
    //mesh.Debug();

    const P1_3d elem;                      // P1 tetrahedron
    FEM_Matrix_2 SK(mesh, elem);           // allocate matrix
    vector<double> uv(SK.Nrows(), 0.0);    //          temperature
    vector<double> fv(SK.Nrows(), 0.0);    //          r.h.s. 
    
    SK.CalculateLaplace(fv,rhs_lap3);      // Fill stiffness matrix and rhs ()
    
    mesh.SetValues(uv, sol_lap3);          // exact solution
    SK.ApplyDirichletBC_Box(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); // brute force Dirichlet via box
    
    SK.Debug();

    vector<double> r_old_n(SK.Nrows(), 0.0);
    vector<double> u_old_n(SK.Nrows(), 0.0);
    vector<double> v_old_n(SK.Nrows(), 0.0);
    vector<double> w_old_n(SK.Nrows(), 0.0);

    double t = 0.5, dt = 0.5;
    //double c = 0.2;
    for (int k = 0; k < t / dt; ++k) {
        vector<double> r_old_m(SK.Nrows(), 0.0);
        vector<double> u_old_m(SK.Nrows(), 0.0);
        vector<double> v_old_m(SK.Nrows(), 0.0);
        vector<double> w_old_m(SK.Nrows(), 0.0);
        for (int n = 0; n <= 3; ++n) {
            //SK1.CalculateLaplace_heat_equation(fv, u_old, dt, t, c);
            SK.CalculateLaplace(fv);
            SK.Debug();


            // Two ways to initialize the vector
            //mesh.SetValues(uv,f_zero);             // user function
            mesh.SetValues(uv, [](double x, double y, double z) -> double {return x *y *z;} ); // lambda function

            auto const save_sol(uv);

            //SK.ApplyDirichletBC(uv,fv);
            SK.ApplyDirichletBC_Box(uv, fv, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
            //SK.Compare2Old(nnode, id, ik, sk);
            //SK.Debug();

            //double tstart = clock();                        // timing
            double tstart = omp_get_wtime();                  // OpenMP

            JacobiSolve(SK, fv, uv );          // solve the system of equations


            //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
            double t1 = omp_get_wtime() - tstart;             // OpenMP
            cout << "JacobiSolve: timing in sec. : " << t1 << endl;


            //mesh.Write_ascii_matlab("uv.txt", uv);
            //mesh.Visualize(uv);

            if (!vectorsAreEqual(uv, save_sol, 1e-3, 10)) {
                cout << "Error in solution." << endl;
            }
        }
        auto s = std::to_string(k);

        mesh.Write_ascii_paraview("uv" + s + ".vtk", uv);
        mesh.Visualize_paraview(uv);
    }

    return true;
}

bool check_P1_2vec_3d()
{     
    Mesh mesh("../salman_problem/NS_3D.txt");
    //Mesh const mesh("../salman_problem/GH_NS_3D.txt");
    mesh.liftToQuadratic();                // needed for P2 elements
    //mesh.Debug();

    const P1_2vec_3d elem;                 // P1-P2vec tetrahedron
    FEM_Matrix_2 SK(mesh, elem);           // allocate matrix
    vector<double> uv(SK.Nrows(), 0.0);    //          temperature
    vector<double> fv(SK.Nrows(), 0.0);    //          r.h.s. 
    
    SK.CalculateLaplace(fv,rhs_lap3);      // Fill stiffness matrix and rhs ()
    
    //mesh.SetValues(uv, sol_lap3);          // exact solution
    //SK.ApplyDirichletBC_Box(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); // brute force Dirichlet via box
    
    SK.Debug();

    //vector<double> r_old_n(SK.Nrows(), 0.0);
    //vector<double> u_old_n(SK.Nrows(), 0.0);
    //vector<double> v_old_n(SK.Nrows(), 0.0);
    //vector<double> w_old_n(SK.Nrows(), 0.0);

    //double t = 0.5, dt = 0.5;
    ////double c = 0.2;
    //for (int k = 0; k < t / dt; ++k) {
        //vector<double> r_old_m(SK.Nrows(), 0.0);
        //vector<double> u_old_m(SK.Nrows(), 0.0);
        //vector<double> v_old_m(SK.Nrows(), 0.0);
        //vector<double> w_old_m(SK.Nrows(), 0.0);
        //for (int n = 0; n <= 3; ++n) {
            ////SK1.CalculateLaplace_heat_equation(fv, u_old, dt, t, c);
            //SK.CalculateLaplace(fv);
            //SK.Debug();


            //// Two ways to initialize the vector
            ////mesh.SetValues(uv,f_zero);             // user function
            //mesh.SetValues(uv, [](double x, double y, double z) -> double {return x *y *z;} ); // lambda function

            //auto const save_sol(uv);

            ////SK.ApplyDirichletBC(uv,fv);
            //SK.ApplyDirichletBC_Box(uv, fv, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
            ////SK.Compare2Old(nnode, id, ik, sk);
            ////SK.Debug();

            ////double tstart = clock();                        // timing
            //double tstart = omp_get_wtime();                  // OpenMP

            //JacobiSolve(SK, fv, uv );          // solve the system of equations


            ////double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
            //double t1 = omp_get_wtime() - tstart;             // OpenMP
            //cout << "JacobiSolve: timing in sec. : " << t1 << endl;


            ////mesh.Write_ascii_matlab("uv.txt", uv);
            ////mesh.Visualize(uv);

            //if (!vectorsAreEqual(uv, save_sol, 1e-3, 10)) {
                //cout << "Error in solution." << endl;
            //}
        //}
        //auto s = std::to_string(k);

        //mesh.Write_ascii_paraview("uv" + s + ".vtk", uv);
        //mesh.Visualize_paraview(uv);
    //}

    return true;
}

