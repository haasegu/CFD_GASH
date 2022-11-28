//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "elements.h"                   // GH: finite elements collection
#include "geom.h"
//#include "geom3.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "userset.h"
#include "vdop.h"
#include "sa_elements.h"
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
bool check_P1_2d();


/**
 * Checks correctness for 3D xxx (in a box domain)
 * using tetrahedrons with linear test function for scalar values
 * and quadratic test functions for the vector values
 *
 * @see P1_2vec_3d
 */




bool check_P1_2vec_2d();

int main(int, char **)
{
    //check_P1_3d();
    check_P1_2vec_2d();

    return 0;
}




/*
bool check_P1_3d()
{
    //Mesh_3d_P1_matlab const mesh("../salman_problem/NS_3D.txt");  // old
    Mesh const mesh("../salman_problem/NS_3D.txt");               // new
    //Mesh const mesh("../salman_problem/GH_NS_3D.txt");
    //mesh.Debug();

    const P1_3d elem;                      // P1 tetrahedron
    //const P1_2vec_3d elem;                      // P2 tetrahedron
    FEM_Matrix_2 SK(mesh, elem);           // allocate matrix
    //int const nP1 = elem.getP1_nodes();
    //int const nP2 = elem.getP2_nodes();
    vector<double> uv(SK.Nrows(), 0.0);    //          temperature
    vector<double> fv(SK.Nrows(), 0.0);    //          r.h.s.


    double t = 0.05, dt = 0.5;
    double c = 0.2;
//    SK.CalculateLaplace(fv,rhs_lap3);      // Fill stiffness matrix and rhs ()
    //SK.CalculateLaplace_heat_equation(fv, uv, dt,t, c);


    mesh.SetValues(uv, sol_lap3);          // exact solution
    SK.ApplyDirichletBC_Box(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); // brute force Dirichlet via box

    SK.Debug();

    vector<double> r_old_n(SK.Nrows(), 0.0);
    vector<double> u_old_n(SK.Nrows(), 0.0);
    vector<double> v_old_n(SK.Nrows(), 0.0);
    vector<double> w_old_n(SK.Nrows(), 0.0);


    for (int k = 0; k < t / dt; ++k) {
        vector<double> r_old_m(SK.Nrows(), 0.0);
        vector<double> u_old_m(SK.Nrows(), 0.0);
        vector<double> v_old_m(SK.Nrows(), 0.0);
        vector<double> w_old_m(SK.Nrows(), 0.0);
        for (int n = 0; n <= 1; ++n) {
            //SK.CalculateLaplace_heat_equation(fv, u_old_m, dt, t, c);
            //SK.CalculateLaplace(fv);
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
*/
/*
bool check_P1_2vec_3d()
{
    Mesh mesh("../salman_problem/NS_3D.txt");
    //Mesh const mesh("../salman_problem/GH_NS_3D.txt");
    mesh.liftToQuadratic();                // needed for P2 elements
    //mesh.Debug();

    const P1_2vec_3d elem;                 // P1-P2vec tetrahedron
    FEM_Matrix_2 SK(mesh, elem);           // allocate matrix
    //int const nP1 = elem.getP1_nodes();
    //int const nP2 = elem.getP2_nodes();
    vector<double> uv(SK.Nrows(), 0.0);    //          temperature
    vector<double> fv(SK.Nrows(), 0.0);    //          r.h.s.

    //std::cout << nP1 << "  " << nP2  << "   " << SK.Nrows() << std::endl;
    vector<double> r(SK.Nrows(), 0.0);
    vector<double> u(SK.Nrows(), 0.0);
    vector<double> v(SK.Nrows(), 0.0);
    vector<double> w(SK.Nrows(), 0.0);
    vector<double> ww(SK.Nrows(), 0.0);

    vector<double> r_ex(SK.Nrows(), 0.0);
    vector<double> u_ex(SK.Nrows(), 0.0);
    vector<double> v_ex(SK.Nrows(), 0.0);
    vector<double> w_ex(SK.Nrows(), 0.0);
    vector<double> uv_ex(SK.Nrows(), 0.0);
    vector<double> uv_bc(SK.Nrows(), 0.0);
    vector<double> uv_ic(SK.Nrows(), 0.0);

    vector<double> r_old_n(SK.Nrows(), 0.0);
    vector<double> u_old_n(SK.Nrows(), 0.0);
    vector<double> v_old_n(SK.Nrows(), 0.0);
    vector<double> w_old_n(SK.Nrows(), 0.0);
    vector<double> r_old_m(SK.Nrows(), 0.0);
    vector<double> u_old_m(SK.Nrows(), 0.0);
    vector<double> v_old_m(SK.Nrows(), 0.0);
    vector<double> w_old_m(SK.Nrows(), 0.0);
    double const t_0 = 0;
    double const dt = 0.001;
    double const t_end = 0.1;
    double const mu = 1.0;
    double const lambda = 1.0;
    double const kp = 0.0;
    double const relax_newton = 0;


    //SK.Navier_Stokes(fv,r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);      // Fill stiffness matrix and rhs ()

    //mesh.SetValues(uv, sol_lap3);          // exact solution
    //SK.ApplyDirichletBC_Box(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); // brute force Dirichlet via box

    //SK.Debug();

    //// Two ways to initialize the vector
    ////mesh.SetValues(uv,f_zero);             // user function


    //mesh.SetValues(uv, [](double x, double y, double z) -> double {return 0.0*x*y;} );  // lambda function

    SK.SetValues(uv_ic, [t_0](double x, double y, double z) -> double {return -(M_PI_4 *M_PI_4 / 2) * (exp(2 * M_PI_4 * x) + exp(2 * M_PI_4 * y) + exp(2 * M_PI_4 * z) + 2 * sin(M_PI_4 * x + M_PI_2 * y) * cos(M_PI_4 * z + M_PI_2 * x) * exp(M_PI_4 * (y + z)) + 2 * sin(M_PI_4 * y + M_PI_2 * z) * cos(M_PI_4 * x + M_PI_2 * y) * exp(M_PI_4 * (z + x)) + 2 * sin(M_PI_4 * z + M_PI_2 * x) * cos(M_PI_4 * y + M_PI_2 * z) * exp(M_PI_4 * (x + y))) * exp(-2 * M_PI_2 *M_PI_2 * t_0);}, // fs     //x*y*z+x*x*x*y*y*y*z-5.0/32
                        [t_0](double x, double y, double z) -> double {return -M_PI_4 * (sin(M_PI_4 * y + M_PI_2 * z) * exp(M_PI_4 * x) + cos(M_PI_4 * x + M_PI_2 * y) * exp(M_PI_4 * z)) * exp(-M_PI_2 *M_PI_2 * t_0);}, // fx   cos(x) *sin(y) *sin(z)    //x+x*x+x*y+x*x*x*y
                        [t_0](double x, double y, double z) -> double {return -M_PI_4 * (sin(M_PI_4 * z + M_PI_2 * x) * exp(M_PI_4 * y) + cos(M_PI_4 * y + M_PI_2 * z) * exp(M_PI_4 * x)) * exp(-M_PI_2 *M_PI_2 * t_0);}, // fy    -sin(x) *cos(y) *sin(z)   //y+x*y+y*y+x*x*y*y
                        [t_0](double x, double y, double z) -> double {return -M_PI_4 * (sin(M_PI_4 * x + M_PI_2 * y) * exp(M_PI_4 * z) + cos(M_PI_4 * z + M_PI_2 * x) * exp(M_PI_4 * y)) * exp(-M_PI_2 *M_PI_2 * t_0);} ); // fz   //-2*z-3*x*z-3*y*z-5*x*x*y*z
   
    SK.Extract_Skalar_Vectors(uv_ic, r, u, v, w);

    for (int ni = 0; ni < t_end / dt; ++ni) {

        double t_ni = t_0 + ni*dt;
        
        SK.SetValues(uv_ex, [t_ni](double x, double y, double z) -> double {return -(M_PI_4 *M_PI_4 / 2) * (exp(2 * M_PI_4 * x) + exp(2 * M_PI_4 * y) + exp(2 * M_PI_4 * z) + 2 * sin(M_PI_4 * x + M_PI_2 * y) * cos(M_PI_4 * z + M_PI_2 * x) * exp(M_PI_4 * (y + z)) + 2 * sin(M_PI_4 * y + M_PI_2 * z) * cos(M_PI_4 * x + M_PI_2 * y) * exp(M_PI_4 * (z + x)) + 2 * sin(M_PI_4 * z + M_PI_2 * x) * cos(M_PI_4 * y + M_PI_2 * z) * exp(M_PI_4 * (x + y))) * exp(-2 * M_PI_2 *M_PI_2 * t_ni);}, // fs     //x*y*z+x*x*x*y*y*y*z-5.0/32
                            [t_ni](double x, double y, double z) -> double {return -M_PI_4 * (sin(M_PI_4 * y + M_PI_2 * z) * exp(M_PI_4 * x) + cos(M_PI_4 * x + M_PI_2 * y) * exp(M_PI_4 * z)) * exp(-M_PI_2 *M_PI_2 * t_ni);}, // fx   cos(x) *sin(y) *sin(z)    //x+x*x+x*y+x*x*x*y
                            [t_ni](double x, double y, double z) -> double {return -M_PI_4 * (sin(M_PI_4 * z + M_PI_2 * x) * exp(M_PI_4 * y) + cos(M_PI_4 * y + M_PI_2 * z) * exp(M_PI_4 * x)) * exp(-M_PI_2 *M_PI_2 * t_ni);}, // fy    -sin(x) *cos(y) *sin(z)   //y+x*y+y*y+x*x*y*y
                            [t_ni](double x, double y, double z) -> double {return -M_PI_4 * (sin(M_PI_4 * x + M_PI_2 * y) * exp(M_PI_4 * z) + cos(M_PI_4 * z + M_PI_2 * x) * exp(M_PI_4 * y)) * exp(-M_PI_2 *M_PI_2 * t_ni);} ); // fz   //-2*z-3*x*z-3*y*z-5*x*x*y*z

            r_old_n = r;
            u_old_n = u;
            v_old_n = v;
            w_old_n = w;
        
        double Error_Newton=1e300;
        for (int mj = 0; mj < 1; ++mj) {
            //do {
            r_old_m = r;
            u_old_m = u;
            v_old_m = v;
            w_old_m = w;
            
            SK.Navier_Stokes(fv, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, mu, lambda, kp);     // Fill stiffness matrix and rhs ()
            //SK.CalculateLaplace(fv);
            //SK.Debug();

            //vector<double> diag(SK.Nrows(), -123456);
            //SK.GetDiag(diag);
            //cout << "diag : " << diag << endl;

            //SK.ApplyDirichletBC(uv,fv);
            // Zero Dirchlet B.C  for update of solution
            SK.ApplyDirichletBC_Box(uv_bc, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            //SK.ApplyDirichletBC_Box(uv_ex, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            //cout << "AA\n";
            //SK.ApplyPeriodicBC_Box_xy(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            //cout << "BB\n";
            //SK.ApplyPeriodicBC_Box_xz(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            //cout << "CC\n";
            //SK.ApplyPeriodicBC_Box_yz(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            cout << "DD\n";

            //SK.Compare2Old(nnode, id, ik, sk);
            //SK.Debug();

            //double tstart = clock();                        // timing
            double tstart = omp_get_wtime();                  // OpenMP

            //JacobiSolve(SK, fv, uv );          // solve the system of equations
            uv = SK.Solve_superlu(fv);
            cout << "EE\n";
            //SK.Mult(ww,uv);
            SK.Defect(ww, fv, uv);
            //cout<<"ww = "<<ww<<endl;
            double sc2 = std::sqrt(L2_scapr(ww, ww));
            cout << "time:  " << t_ni << "   ||defect of direct solver||  " << sc2 << endl;

            //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
            double t1 = omp_get_wtime() - tstart;             // OpenMP
            cout << "Solve: timing in sec. : " << t1 << endl;

            ////mesh.Write_ascii_matlab("uv.txt", uv);
            ////mesh.Visualize(uv);

            if (!vectorsAreEqual(uv, uv_ex, 1e-5, 10)) {
                cout << "Error in solution." << endl;
            }

            SK.Extract_Skalar_Vectors(uv, r, u, v, w);
            
            

            for (size_t i = 0; i < r.size(); ++i) {
                r[i] += relax_newton*r_old_m[i];
            }
            for (size_t j = 0; j < u.size(); ++j) {
                u[j] += relax_newton*u_old_m[j];
                v[j] += relax_newton*v_old_m[j];
                w[j] += relax_newton*w_old_m[j];
            }

            cout << " uv: " << uv << endl;
            cout<<"ccccccccccccccccccccccccccccccccccc"<<endl;
            cout << " uv_ex: " << uv_ex << endl;

            Error_Newton = SK.Error_Newton_Step(r_old_m, u_old_m, v_old_m, w_old_m, r, u, v, w);
            cout << "%%%%%%%%%%%%  Error_Newton_Step = " << Error_Newton << endl;

        }//while(Error_Newton>1e1);

        //ww = getAbsError(uv, uv_ex);
       //cout << "Error_EX_Nu = " << ww << endl;

        //CompareVectors(r, e_node, r_ex, 1e-6);

        //vectorsAreEqual(u, u_ex, 1e-6, 10)<<endl;

        //JacobiSolve(SK, fv, uv );          // solve the system of equations
        //SK.Extract_Skalar_Vectors(uv, r, u, v, w);

        //auto s = std::to_string(n);

        //SK.Write_ascii_paraview_r("w_r" + s + ".vtk", r);
        //SK.Visualize_paraview_r(r);
        //SK.Write_ascii_paraview_u("w_u" + s + ".vtk", u);
        //SK.Visualize_paraview_u(u);
        //SK.Write_ascii_paraview_v("w_v" + s + ".vtk", v);
        //SK.Visualize_paraview_v(v);
        //SK.Write_ascii_paraview_w("w_w" + s + ".vtk", w);
        //SK.Visualize_paraview_w(w);
    }

    //cout << "nP1 = "<<nP1<<", "<<"nP2 = "<<nP2<<endl;
    //cout << "r_old_n: " << r_old_n << endl;
    //cout << "r_old_m: " << r_old_m << endl;
    //cout << "u_old_n: " << u_old_n << endl;
    //cout << "u_old_m: " << u_old_m << endl;
    //cout << "v_old_n: " << v_old_n << endl;
    //cout << "v_old_m: " << v_old_m << endl;
    //cout << "w_old_n: " << w_old_n << endl;
    //cout << "w_old_m: " << w_old_m << endl;

    //SK.Extract_Skalar_Vectors(uv, r, u, v, w);

    //auto s = std::to_string(n);
    //SK.Write_ascii_paraview_r("w_r" + s + ".vtk", r);

    //cout << "uv: " << uv << endl;
    cout << " r: " << r << endl;
    cout << "u: " << u << endl;
    cout << "v: " << v << endl;
    cout << "w: " << w << endl;
    SK.Visualize_paraview_r(r);

    //SK.Write_ascii_paraview_u("w_u" + s + ".vtk", u);
    //SK.Visualize_paraview_u(u);
    //SK.Write_ascii_paraview_v("w_v" + s + ".vtk", v);
    //SK.Visualize_paraview_v(v);
    //SK.Write_ascii_paraview_w("w_w" + s + ".vtk", w);
    //SK.Visualize_paraview_w(w);

    //double tstart = clock();                        // timing
    double tstart = omp_get_wtime();                  // OpenMP

    //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
    double t1 = omp_get_wtime() - tstart;             // OpenMP
    cout << "JacobiSolve: timing in sec. : " << t1 << endl;
*/
bool check_P1_2vec_2d()
{
    Mesh mesh("square.txt");
    //Mesh const mesh("../salman_problem/GH_NS_3D.txt");
    mesh.liftToQuadratic();                // needed for P2 elements
    //mesh.Debug();

    const P1_2vec_2d elem;                 // P1-P2vec tetrahedron
    FEM_Matrix_2<P1_2vec_2d> SK(mesh, elem);           // allocate matrix
    //int const nP1 = elem.getP1_nodes();
    //int const nP2 = elem.getP2_nodes();
    vector<double> uv(SK.Nrows(), 0.0);    //          temperature
    vector<double> fv(SK.Nrows(), 0.0);    //          r.h.s.

    //std::cout << nP1 << "  " << nP2  << "   " << SK.Nrows() << std::endl;
    vector<double> r(SK.Nrows(), 0.0);
    vector<double> u(SK.Nrows(), 0.0);
    vector<double> v(SK.Nrows(), 0.0);
    vector<double> w(SK.Nrows(), 0.0);
    vector<double> ww(SK.Nrows(), 0.0);

    vector<double> r_ex(SK.Nrows(), 0.0);
    vector<double> u_ex(SK.Nrows(), 0.0);
    vector<double> v_ex(SK.Nrows(), 0.0);
    vector<double> w_ex(SK.Nrows(), 0.0);
    
    vector<double> r_r(SK.Nrows(), 0.0);
    vector<double> u_u(SK.Nrows(), 0.0);
    vector<double> v_v(SK.Nrows(), 0.0);
    vector<double> w_w(SK.Nrows(), 0.0);
    
    vector<double> uv_ex(SK.Nrows(), 0.0);
    vector<double> uv_bc(SK.Nrows(), 0.0);
    vector<double> uv_ic(SK.Nrows(), 0.0);
    vector<double> uv_out(SK.Nrows(), 0.0);
    vector<double> fv_ex(SK.Nrows(), 0.0);


    vector<double> r_old_n(SK.Nrows(), 0.0);
    vector<double> u_old_n(SK.Nrows(), 0.0);
    vector<double> v_old_n(SK.Nrows(), 0.0);
    vector<double> w_old_n(SK.Nrows(), 0.0);
    vector<double> r_old_m(SK.Nrows(), 0.0);
    vector<double> u_old_m(SK.Nrows(), 0.0);
    vector<double> v_old_m(SK.Nrows(), 0.0);
    vector<double> w_old_m(SK.Nrows(), 0.0);
    double const t_0 = 0;
    double const dt = 0.0001;
    double const t_end = 1;
    double const mu = 0.1;
    double const lambda = 0.1;
    double const r0 = 0.2;
    double const a = 1.0;
    double const gamma = 1.0;
    double const kt = 0.1;
    double const cp = 0.1;
    double const relax_newton = 1.0;
    
    

       SK.SetValues(uv_ic,  [t_0](double x, double y) -> double {return  1-0.5*tanh(y-0.5);},//1;}, 
                            [t_0](double x, double y) -> double {return  sin(M_PI*x)*sin(M_PI*x)*sin(2*M_PI*y);}, //-cos(M_PI*x)*sin(M_PI*y)*exp(-2*M_PI*t_0);},
                            [t_0](double x, double y) -> double {return  -sin(M_PI*y)*sin(M_PI*y)*sin(2*M_PI*x);},
                            [t_0,r0](double x, double y) -> double {  if(0<=sqrt(pow(x-0.5,2)+pow(y-0.5,2))<r0/2){ return 2.0*sqrt(pow(x-0.5,2)+pow(y-0.5,2))/r0;} else if(r0/2<=sqrt(pow(x-0.5,2)+pow(y-0.5,2))<r0)
								                     {return 2.0*(1-sqrt(pow(x-0.5,2)+pow(y-0.5,2))/r0);} else if(r0<=sqrt(pow(x-0.5,2)+pow(y-0.5,2))<=0.5){return 0;} return 0;}); //sin(M_PI*x)*cos(M_PI*y)*exp(-2*M_PI*t_0);});

                            
        SK.Extract_Skalar_Vectors(uv_ic, r, u, v, w);

        for (int ni = 0; ni < t_end / dt; ++ni) {

        double t_ni = t_0 + ni*dt;
        
            r_old_n = r;
            u_old_n = u;
            v_old_n = v;
            w_old_n = w;
                            
        double Error_Newton=1e300;
        for (int mj = 0; mj < 5; ++mj) {
            //do {
            r_old_m = r;
            u_old_m = u;
            v_old_m = v;
            w_old_m = w;
            
            
            SK.Navier_Stokes(fv, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, mu, lambda, a, gamma, t_ni, kt, cp);     // Fill stiffness matrix and rhs ()
            //cout << "fv= " <<fv << endl;
            //SK.CalculateLaplace(fv);
            //SK.Debug();

            //vector<double> diag(SK.Nrows(), -123456);
            //SK.GetDiag(diag);
            //cout << "diag : " << diag << endl;

            //SK.ApplyDirichletBC(uv,fv);
            // Zero Dirchlet B.C  for update of solution
            //SK.ApplyDirichletBC_Box(uv, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            SK.ApplyDirichletBC_Box(uv_out, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0);
            //cout << "AA\n";
            //SK.ApplyPeriodicBC_Box_xy(uv_ex, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            //cout << "BB\n";
            //SK.ApplyPeriodicBC_Box_xz(uv_ex, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            //cout << "CC\n";
            //SK.ApplyPeriodicBC_Box_yz(uv_ex, fv, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            cout << "DD\n";

            //SK.Compare2Old(nnode, id, ik, sk);
            //SK.Debug();

            //double tstart = clock();                        // timing
            double tstart = omp_get_wtime();                  // OpenMP

            //JacobiSolve(SK, fv, uv );          // solve the system of equations
            uv = SK.Solve_superlu(fv);
            cout << "EE\n";
            //SK.Mult(ww,uv);
            //SK.Defect(ww, fv, uv);
            //cout<<"ww = "<<ww<<endl;
            double sc2 = std::sqrt(L2_scapr(ww, ww));
            cout << "time:  " << t_ni << "   ||defect of direct solver||  " << sc2 << endl;

            //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
            double t1 = omp_get_wtime() - tstart;             // OpenMP
            cout << "Solve: timing in sec. : " << t1 << endl;

            if (!vectorsAreEqual(uv, uv_ex, 1e-5, 10)) {
                cout << "Error in solution." << endl;
            }

            SK.Extract_Skalar_Vectors(uv_out, r, u, v, w);
            
            SK.Extract_Skalar_Vectors(uv, r_r, u_u, v_v, w_w);
            

            for (size_t i = 0; i < r.size(); ++i) {
                r[i] = r_old_m[i]+relax_newton*r_r[i];
            }
            for (size_t j = 0; j < u.size(); ++j) {
                u[j] = u_old_m[j]+relax_newton*u_u[j];
                v[j] = v_old_m[j]+relax_newton*v_v[j];
                w[j] = w_old_m[j]+relax_newton*w_w[j];
            }
            
            
            //SK.Extract_Skalar_Vectors(uv_ex, r_r, u_u, v_v, w_w);
            
            //ww = getAbsError(u,u_u);
            //cout << "Error_EX_Nu = " << ww << endl;
            
            
            //cout << "uv: " << uv << endl;
            //cout << "r: " << r << endl;
            //cout << "u: " << u << endl;
            //cout << "v: " << v << endl;
            
            //auto s = std::to_string(n);
            //SK.Write_ascii_paraview_l("w_r.vtk", r);
            //SK.Visualize_paraview_l(r);
            
            

            

        }//while(Error_Newton>1e1);

       //ww = getAbsError(uv, uv_ex);
       //cout << "Error_EX_Nu = " << ww << endl;

        //CompareVectors(r, e_node, r_ex, 1e-6);

        //vectorsAreEqual(u, u_ex, 1e-6, 10)<<endl;

        //JacobiSolve(SK, fv, uv );          // solve the system of equations
        //SK.Extract_Skalar_Vectors(uv, r, u, v, w);
        
        
        //cout << "uv: " << uv << endl;
        cout << "r: " << r << endl;
        //cout << "u: " << u << endl;
        //cout << "v: " << v << endl;
        //cout << "w: " << w << endl;
        //SK.Write_ascii_paraview_l_r("w_r.vtk", r);
        //SK.Visualize_paraview_l_r(r);
        
        //SK.Write_ascii_paraview_q_w("w_w.vtk", w);
        //SK.Visualize_paraview_q_w(w);
        
        auto s = std::to_string(ni);
        if(0.0==t_ni ||0.01==t_ni ||0.05==t_ni ||0.1==t_ni ||0.5==t_ni ||1.0==t_ni){
          SK.Write_ascii_paraview_l_r("w_r_"+s+".vtk", r);
          //SK.Visualize_paraview_l_r(r);
          
          SK.Write_ascii_paraview_q_u("w_u_"+s+".vtk", u);
          SK.Visualize_paraview_q_u(u);
          
          SK.Write_ascii_paraview_q_v("w_v_"+s+".vtk", v);
          //SK.Visualize_paraview_q_v(v);
          
          SK.Write_ascii_paraview_q_w("w_w_"+s+".vtk", w);
          //SK.Visualize_paraview_q_w(w);
	  }

        
    }
    
          SK.Write_ascii_paraview_l_r("w_r.vtk", r);
          //SK.Visualize_paraview_l_r(r);
          
          SK.Write_ascii_paraview_q_u("w_u.vtk", u);
          //SK.Visualize_paraview_q_u(u);
          
          SK.Write_ascii_paraview_q_v("w_v.vtk", v);
          //SK.Visualize_paraview_q_v(v);
          
          SK.Write_ascii_paraview_q_w("w_w.vtk", w);
          //SK.Visualize_paraview_q_w(w);

     //cout << "r: " << r << endl;
     //cout << "u: " << u << endl;
     //cout << "v: " << v << endl;

    //double tstart = clock();                        // timing
    double tstart = omp_get_wtime();                  // OpenMP

    //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
    double t1 = omp_get_wtime() - tstart;             // OpenMP
    cout << "JacobiSolve: timing in sec. : " << t1 << endl;

    return true;
}











//########################################################################
/*
    int nthreads;                                  // OpenMP
    #pragma omp parallel default(none) shared(cout,nthreads)
    {
        int const th_id  = omp_get_thread_num();   // OpenMP
        int const nthrds = omp_get_num_threads();  // OpenMP
        stringstream ss;
        ss << "C++: Hello World from thread " << th_id << " / " << nthrds << endl;
        #pragma omp critical
        {
            cout << ss.str();                      // output to a shared ressource
        }
        #pragma omp master
        nthreads = nthrds;                         // transfer nn to to master thread
    }
    cout << "   " << nthreads << "   threads have been started." << endl;
*/
// ##################### STL ###########################################
/*
{
    //Mesh_2d_3_matlab const mesh("square_tiny.txt");
    Mesh_2d_3_matlab const mesh("square_100.txt");
    //Mesh_2d_3_matlab const mesh("L_shape.txt");
    //mesh.Debug();

    CRS_Matrix SK(mesh);                   // CRS matrix
    //SK.Debug();

    vector<double> uv(SK.Nrows(),0.0);     // temperature
    vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

    SK.CalculateLaplace(fv);
    //SK.Debug();

    // Two ways to initialize the vector
    //mesh.SetValues(uv,f_zero);             // user function
    mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function

    SK.ApplyDirichletBC(uv,fv);
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
    }
*/
// ##################### STL ###########################################
/*
{
    //Mesh_2d_3_matlab const mesh("square_tiny.txt");
    Mesh_3d_4_matlab const mesh("../salman_problem/NS_3D.txt");
    //Mesh_2d_3_matlab const mesh("L_shape.txt");
    //mesh.Debug();

     //CRS_Matrix SK(mesh);                   // CRS matrix
     //SK.Debug();

    FEM_Matrix SK1(mesh);
    cout << "########################################################\n";
    SK1.Debug();


    vector<double> uv(SK1.Nrows(),0.0);     // temperature
    vector<double> fv(SK1.Nrows(),0.0);     // r.h.s.
    vector<double> r_old_n(SK1.Nrows(),0.0);
    vector<double> u_old_n(SK1.Nrows(),0.0);
    vector<double> v_old_n(SK1.Nrows(),0.0);
    vector<double> w_old_n(SK1.Nrows(),0.0);




    double t=0.5, dt = 0.5, mu = 0.1, lambda = 0.1, kp = 0.5;
    for (int n = 0; n <= t/dt; ++n)
    {
	vector<double> r_old_m(SK1.Nrows(),0.0);
    vector<double> u_old_m(SK1.Nrows(),0.0);
    vector<double> v_old_m(SK1.Nrows(),0.0);
    vector<double> w_old_m(SK1.Nrows(),0.0);

    for(int m=0; m<=3; ++m){
    //SK1.CalculateLaplace_heat_equation(fv, u_old, dt, t, c);
    SK1.CalculateLaplace(fv);
    SK1.Debug();


    // Two ways to initialize the vector
    //mesh.SetValues(uv,f_zero);             // user function
    mesh.SetValues(uv, [](double x, double y, double z) -> double {return x*y*z;} );  // lambda function

    auto const save_sol(uv);

    //SK1.ApplyDirichletBC(uv,fv);
    SK1.ApplyDirichletBC_Box(uv,fv,1.0,1.0,1.0,1.0,1.0,1.0);
    //SK.Compare2Old(nnode, id, ik, sk);
    //SK.Debug();

    //double tstart = clock();                        // timing
    double tstart = omp_get_wtime();                  // OpenMP

    JacobiSolve(SK1, fv, uv );          // solve the system of equations






    //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
    double t1 = omp_get_wtime() - tstart;             // OpenMP
    cout << "JacobiSolve: timing in sec. : " << t1 << endl;


    //mesh.Write_ascii_matlab("uv.txt", uv);
    //mesh.Visualize(uv);

    if (!vectorsAreEqual(uv,save_sol,1e-3,10))
    {
		cout << "Error in solution." << endl;
	}
}
	auto s = std::to_string(n);

    mesh.Write_ascii_paraview("uv"+s+".vtk", uv);
    mesh.Visualize_paraview(uv);
    }

}
*/
/*if(0==ni){
			   for(unsigned int i=0; i<r.size(); ++i)
                  {
                    r_old_n[i] = uv[i];
                  }
               for(unsigned int j=0; j<u.size(); ++j)
                  {
					int jk = r.size()+3*j;
                    u_old_n[j] = uv[jk+0];
                    v_old_n[j] = uv[jk+1];
                    w_old_n[j] = uv[jk+2];
                  }

			   	}
		    else{
			   for(unsigned int i=0; i<r.size(); ++i)
                  {
                    r_old_n[i] = r[i];
                  }
               for(unsigned int j=0; j<u.size(); ++j)
                  {
                    u_old_n[j] = u[j];
                    v_old_n[j] = v[j];
                    w_old_n[j] = w[j];
                  }

			   	}*/

/*if(0==ni && 0==mj){
               for(unsigned int i=0; i<r.size(); ++i)
                  {
                    r_old_m[i] = uv[i];
                  }
               for(unsigned int j=0; j<u.size(); ++j)
                  {
                    u_old_m[j] = uv[r.size()+j];
                    v_old_m[j] = uv[r.size()+u.size()+j];
                    w_old_m[j] = uv[r.size()+2*u.size()+j];
                  }

               	}
if(0==mj){
       for(unsigned int i=0; i<r.size(); ++i)
          {
            r_old_m[i] = r_old_n[i];
          }
       for(unsigned int j=0; j<u.size(); ++j)
          {
            u_old_m[j] = u_old_n[j];
            v_old_m[j] = v_old_n[j];
            w_old_m[j] = w_old_n[j];
          }

       	}
else{
   for(unsigned int i=0; i<r.size(); ++i)
      {
        r_old_m[i] = r[i];
      }
   for(unsigned int j=0; j<u.size(); ++j)
      {
        u_old_m[j] = u[j];
        v_old_m[j] = v[j];
        w_old_m[j] = w[j];
      }

   	}*/
