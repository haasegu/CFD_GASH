//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "geom.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "userset.h"
#include "vdop.h"

#include "omp.h"
#include <chrono>           // timing
#include <cmath>
#include <iostream>
using namespace std;
using namespace std::chrono;  // timing


void Test_solver(int lev, vector<int> const & nthreads, Multigrid& ggm);

int main(int argc , char **argv )
{
    const int numprocs = 1;
    const int myrank   = 0;

    if (myrank == 0)
    {
        cout << "\n There are " << numprocs << " processes running.\n \n";
    }

    const auto procx = static_cast<int>(sqrt(numprocs + 0.0));
    const int  procy = procx;

    if (procy * procx != numprocs)
    {
        cout << "\n Wrong number of processors !\n \n";
        return -1;
    }
    else
    {
// #####################################################################
//      Here starts the real code
// #####################################################################
        ////bool ScaleUp = !true;
        //int nx, ny, NXglob, NYglob; /* number of local intervals on (xl,xr)=:nx, (yb,yt)=:ny */
        ////nx = 1024;
        ////ny = 1024;
        //nx = 100;
        //ny = 100;
        //NXglob = nx * procx;
        //NYglob = ny * procy;
        //cout << "Intervalls: " << NXglob << " x " << NYglob << endl;

// ##################### STL ###########################################
//{
        //Mesh_2d_3_square const mesh(nx, ny);
        ////mesh.Debug();

        //FEM_Matrix SK(mesh);                   // CRS matrix
        ////SK.Debug();
        //vector<double> uv(SK.Nrows(),0.0);     // temperature
        //vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

        //SK.CalculateLaplace(fv);
        ////SK.Debug();

        ////mesh.SetU(uv);         // deprecated
        ////mesh.SetF(fv);         // deprecated
        //// Two ways to initialize the vector
        ////mesh.SetValues(uv,f_zero);             // functional
        //mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function

        //SK.ApplyDirichletBC(uv,fv);
        ////SK.Compare2Old(nnode, id, ik, sk);
        ////SK.Debug();

        //double tstart = clock();                        // timing

        //JacobiSolve(SK, fv, uv );          // solve the system of equations

        //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
        //cout << "JacobiSolve: timing in sec. : " << t1 << endl;

        ////CompareVectors(uv, nnode, u, 1e-6);    // Check correctness

        ////mesh.SaveVectorP("t.dat", uv);
        ////mesh.Visualize(uv);

//}
// ##################### STL ###########################################
#define MG
#ifndef MG
{
        int nrefine = 3;
        if (argc>1)  nrefine = atoi(argv[1]);
        //Mesh const mesh("square_tiny.txt");
        //Mesh const mesh("square_100.txt");
        //Mesh const mesh("L_shape.txt");

        //Mesh const mesh_c("square_tiny.txt");
        Mesh const mesh_c("../salman_problem/NS_3D.txt");
        //Mesh const mesh_c("square.txt");
        mesh_c.checkObtuseAngles();
        //return 0;       
        //mesh_c.Debug();
        //mesh_c.DebugEdgeBased();

        //RefinedMesh mesh(mesh_c);                  // OK, works
        //Mesh const mesh("square_tiny.txt");
        ////mesh.Debug();
        //mesh.RefineAllElements(nrefine);           // OK, works

        gMesh_Hierarchy ggm(mesh_c,nrefine);
        const Mesh& mesh=ggm.finest();

        //mesh.Debug();
        //mesh.DebugEdgeBased();


        FEM_Matrix SK(mesh);                   // CRS matrix
        //SK.writeBinary("sparseMatrix.bin");
        //SK.Debug();

        vector<double> uv(SK.Nrows(),0.0);     // temperature
        vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

        SK.CalculateLaplace(fv);
        //SK.CheckRowSum();
        SK.CheckMatrix();
        //return 0;
        //SK.Debug();

        //mesh.SetU(uv);         // deprecated
        // Two ways to initialize the vector
        //mesh.SetValues(uv,f_zero);             // user function
        //mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function
        //mesh.SetValues(uv, [](double x, double y) -> double {return 5e-3*(x+1)*(y+1);} );  // lambda function
        //
        mesh.SetValues(uv, [](double x, double y) -> double {
            return x * x * std::sin(2.5 * M_PI * y);
            } );

        SK.ApplyDirichletBC(uv,fv);
        
        //SK.Compare2Old(nnode, id, ik, sk);
        //SK.Debug();
        auto exact_sol(uv);
        
        //SK.Mult(fv,uv);

        auto t3 = system_clock::now(); // start timer

        JacobiSolve(SK, fv, uv );          // solve the system of equations

        auto t4 = system_clock::now();  // stop timer
        auto duration = duration_cast<microseconds>(t4 - t3);        // duration in microseconds
        double t_diff = static_cast<double>(duration.count()) / 1e6; // overall duration in seconds
        cout << "JacobiSolve: timing in sec. : " << t_diff << endl;
        
        auto [val,idx] = findLargestAbsError(exact_sol, uv, 1e+6, 100);
        
        //mesh.Visualize(getAbsError(exact_sol, uv));
        

        //mesh.Write_ascii_matlab("uv.txt", uv);
        //mesh.Visualize(uv);
    }
#else
    {
        int nrefine = 3;
        if (argc>1)  nrefine = atoi(argv[1]);

        //Multigrid ggm(Mesh("square_tiny.txt"),nrefine);
        Multigrid ggm(Mesh("../salman_problem/NS_3D.txt"),nrefine);

        ggm.DefineOperators();

        cout << "\n#############  SOLVE   ###############\n";

        double tstart = omp_get_wtime();                  // OpenMP

        //ggm.JacobiSolve(my_level);
        //ggm.MG_Step(my_level, 1, true, 1);
        ggm.MG_Solve(2, 1e-6, 1);

        double t1 = omp_get_wtime() - tstart;             // OpenMP
        cout << "MgSolve: timing in sec. : " << t1 << "   for " << ggm.Ndofs()<< " dofs"<< endl;
        

        Test_solver(nrefine-1, {1,2,4,8,16,32,64,128,256}, ggm);

        //int my_level=nrefine-1;
        //const auto &ml=ggm.GetMesh(my_level);
        //const auto &sl=ggm.GetSolution(my_level);
        //ml.Visualize(sl);
        //////ml.Visualize_paraview(sl);
        ////ml.Export_scicomp("level_"+to_string(my_level));
        
        //int my_level=nrefine-1;
        //const auto &mesh=ggm.GetMesh(my_level);
        //const auto &uv=ggm.GetSolution(my_level);
        //vector<double> exact_sol(size(uv)); 
        //mesh.SetValues(exact_sol, [](double x, double y) -> double {
            //return x * x * std::sin(2.5 * M_PI * y);
            //} );
        //mesh.Visualize(getAbsError(exact_sol, uv));
    }
#endif
    return 0;
}
}


void Test_solver(int lev, vector<int> const & nthreads, Multigrid& ggm)
{
    cout << endl << endl << "-------------------------------------" << endl;
    cout << "MgSolve: timing in sec. for " << ggm.Ndofs()<< " dofs"<< endl;
    cout << "sec         threads" << endl;
    vector<double> mg_time(size(nthreads),-1.0);
    
    for (size_t k=0; k<size(nthreads); ++k)
    {
        omp_set_num_threads(nthreads.at(k));
        double tstart = omp_get_wtime();
        ggm.MG_Solve(2, 1e-6, 1);
        double t1 = omp_get_wtime() - tstart;
        mg_time.at(k) = t1;
        cout << t1 << " : " << nthreads.at(k) << endl;
    }
}
