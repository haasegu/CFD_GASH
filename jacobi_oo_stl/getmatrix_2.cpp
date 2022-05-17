#include "getmatrix_2.h"

#include <vector>

using namespace std;

FEM_Matrix_P2::FEM_Matrix_P2(Mesh const &mesh, int ndof_v)
    : FEM_Matrix(mesh,1)
{
    assert(mesh.NverticesElements()==10);          // quadratic tetrahedral elements?
    vector<int> const &ia_geom = mesh.GetConnectivity();   // geometric connectivity
    
    //Derive_Matrix_Pattern();
    //Skalar2VectorMatrix(ndof_v);
    //return;
}


FEM_Matrix_P2::~FEM_Matrix_P2()
{}

