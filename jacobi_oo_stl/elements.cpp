#include "elements.h"
#include "geom.h"
#include "getmatrix.h"
#include "userset.h"
#include "utils.h"
#include "vdop.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;

vector<int> shrinkColumns(vector<int> const &v2, int col2, int col1)
{
    assert( v2.size()%col2 == 0 );
    assert( col1<=col2 );
    int const nrows = static_cast<int>(v2.size()/col2);
    vector<int> v1(nrows*col1);
    
    for (int k=0; k<nrows; ++k)
    {
        for (int d=0; d<col1; ++d)
        {
            v1.at(k*col1+d) = v2.at(k*col2+d);
        }
    }
    
    return v1;
}

vector<int> mergeColumns(vector<int> const &v1, int col1, vector<int> const &v2, int col2)
{
    assert( v2.size()%col2 == 0 );
    assert( v1.size()%col1 == 0 );
    assert( v1.size()/col1 == v2.size()/col2 );
    int const nrows = static_cast<int>(v2.size()/col2);
    int const ncols = col1+col2;
    vector<int> vm(nrows*ncols);
    
    
    for (int k=0; k<nrows; ++k)
    {
        for (int d=0; d<col1; ++d)
        {
            vm.at(k*ncols+d) = v1.at(k*col1+d);       // copy data v1
        }
        for (int d=0; d<col2; ++d)
        {
            vm.at(k*ncols+d+col1) = v2.at(k*col2+d);  // copy data v2
        }
    }
    
    return vm;
}



vector<int> indexVectorBlowUp(vector<int> const &vin, int ndof_v, int offset)
{
    vector<int> wout(vin.size()*ndof_v);
    
    for (size_t k=0; k<vin.size(); ++k)
    {
        int const kv=k*ndof_v;
        for (int d=0; d<ndof_v; ++d)
        {
            wout[kv+d] = ndof_v*vin[k]+offset+d;
        }
    }
    
    return wout;
}

//######################################################################

Element::~Element() {}

vector<int> P1_2d::getGlobalDOFIndices(vector<int> const &ia_geom)  const
{
    assert(0 == nDOFs_loc() % nVertices_loc());
    int const ndof_v = nDOFs_loc()/nVertices_loc();  // DOFs per vertex
    return indexVectorBlowUp(ia_geom,ndof_v,0);
}

vector<int> P1_3d::getGlobalDOFIndices(vector<int> const &ia_geom)  const
{
    assert(0 == nDOFs_loc() % nVertices_loc());
    int const ndof_v = nDOFs_loc()/nVertices_loc();  // DOFs per vertex
    return indexVectorBlowUp(ia_geom,ndof_v,0);
}

vector<int> P2_3d::getGlobalDOFIndices(vector<int> const &ia_geom)  const 
{
    assert(0 == nDOFs_loc() % nVertices_loc());
    int const ndof_v = nDOFs_loc()/nVertices_loc();  // DOFs per vertex
    //cout << "LLLLLLLLLLL   " << ndof_v  << endl;
    return indexVectorBlowUp(ia_geom,ndof_v,0);    
}

vector<int> P1_2vec_3d::getGlobalDOFIndices(vector<int> const &ia_geom_P2)  const
{
    int const nv1(4);              // #vertices tet linear
    int const nv2(10);             // #vertices tet quadratic    
    // first: P1 part - scalar
    vector<int> ia_geom_P1 = shrinkColumns(ia_geom_P2, nv2, nv1);
    auto const ia_p1(indexVectorBlowUp(ia_geom_P1,1,0));
    int const nP1_dofs = *max_element(cbegin(ia_p1),cend(ia_p1));  // offset for next numbering
    
    // second: P2 part - vector
    auto const ia_p2(indexVectorBlowUp(ia_geom_P2,3,nP1_dofs));
    
    // merging of  ia_p1 and ia_p2
    int const col1=nv1*1;
    int const col2=nv2*3;
    vector<int> ia_dofs(mergeColumns(ia_p1, col1, ia_p2, col2));
    
    return ia_dofs;
}

////######################################################################

//  general routine for lin. triangular elements
void P1_2d::CalcLaplace(
    int const ial[3], double const xc[], 
    vector<vector<double>> &ske, vector<double> &fe,
    const function<double(double, double)> &f_func  ) const
{
    assert(nVertices_loc()==3);
    assert(nDOFs_loc()==3);
    //cout << ske.size() << "  " << nDOFs_loc() << endl;
    assert(static_cast<int>(ske.size())==nDOFs_loc());
    const int  i1 = 2*ial[0], i2 = 2*ial[1], i3 = 2*ial[2];
    const double x13 = xc[i3 + 0] - xc[i1 + 0],  y13 = xc[i3 + 1] - xc[i1 + 1],
                 x21 = xc[i1 + 0] - xc[i2 + 0],  y21 = xc[i1 + 1] - xc[i2 + 1],
                 x32 = xc[i2 + 0] - xc[i3 + 0],  y32 = xc[i2 + 1] - xc[i3 + 1];
    const double jac = std::abs(x21 * y13 - x13 * y21);

    ske[0][0] = 0.5 / jac * (y32 * y32 + x32 * x32);
    ske[0][1] = 0.5 / jac * (y13 * y32 + x13 * x32);
    ske[0][2] = 0.5 / jac * (y21 * y32 + x21 * x32);
    ske[1][0] = ske[0][1];
    ske[1][1] = 0.5 / jac * (y13 * y13 + x13 * x13);
    ske[1][2] = 0.5 / jac * (y21 * y13 + x21 * x13);
    ske[2][0] = ske[0][2];
    ske[2][1] = ske[1][2];
    ske[2][2] = 0.5 / jac * (y21 * y21 + x21 * x21);

    const double xm    = (xc[i1 + 0] + xc[i2 + 0] + xc[i3 + 0]) / 3.0,
                 ym    = (xc[i1 + 1] + xc[i2 + 1] + xc[i3 + 1]) / 3.0;
    fe[0] = fe[1] = fe[2] = 0.5 * jac * f_func(xm, ym) / 3.0;
}

void P1_2d::CalcLaplace(
    int const ial[3], double const xc[], 
    vector<vector<double>> &ske, vector<double> &fe) const
{
    CalcLaplace(ial,xc,ske,fe,rhs_lap2);
}

void P1_3d::CalcLaplace(
    int const ial[4], double const xc[], 
    vector<vector<double>> &ske, vector<double> &fe,
    const function<double(double, double, double)> &f_func) const
{
    assert(nVertices_loc()==4);
    assert(nDOFs_loc()==4);
    assert(static_cast<int>(ske.size())==nDOFs_loc());
    assert(3==nDim_loc());
    // global (geom.) vertex indices
    const int  i1 = 3*ial[0], i2 = 3*ial[1], i3 = 3*ial[2], i4 = 3*ial[3]; 
    // coordinates of geom. vertices
    const double x1{xc[i1+0]}, x2{xc[i2+0]}, x3{xc[i3+0]}, x4{xc[i4+0]};
    const double y1{xc[i1+1]}, y2{xc[i2+1]}, y3{xc[i3+1]}, y4{xc[i4+1]};
    const double z1{xc[i1+2]}, z2{xc[i2+2]}, z3{xc[i3+2]}, z4{xc[i4+2]};
    //  center of gravity of the tetrahedron
    const double xm{(x1+x2+x3+x4)/4}, ym{(y1+y2+y3+y4)/4}, zm{(z1+z2+z3+z4)/4};
    
    // copy-paste from Matlab file tet_elem.m
    const double detA=x1*y3*z2 - x1*y2*z3 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1 - x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1 + x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2;

    const double gradPhi[4][3] = { //       1/detA*
        {y3*z2 - y2*z3 + y2*z4 - y4*z2 - y3*z4 + y4*z3, x2*z3 - x3*z2 - x2*z4 + x4*z2 + x3*z4 - x4*z3, x3*y2 - x2*y3 + x2*y4 - x4*y2 - x3*y4 + x4*y3},
        {y1*z3 - y3*z1 - y1*z4 + y4*z1 + y3*z4 - y4*z3, x3*z1 - x1*z3 + x1*z4 - x4*z1 - x3*z4 + x4*z3, x1*y3 - x3*y1 - x1*y4 + x4*y1 + x3*y4 - x4*y3},
        {y2*z1 - y1*z2 + y1*z4 - y4*z1 - y2*z4 + y4*z2, x1*z2 - x2*z1 - x1*z4 + x4*z1 + x2*z4 - x4*z2, x2*y1 - x1*y2 + x1*y4 - x4*y1 - x2*y4 + x4*y2},
        {y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2, x2*z1 - x1*z2 + x1*z3 - x3*z1 - x2*z3 + x3*z2, x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2},
    }; 
    
    for (int row=0; row<nVertices_loc(); ++row)
    {
        for (int col=0; col<nVertices_loc(); ++col)
        {
            ske[row][col] = gradPhi[row][0]*gradPhi[col][0]+gradPhi[row][1]*gradPhi[col][1]+gradPhi[row][2]*gradPhi[col][2];
            ske[row][col] /= (detA*6.0);        // GH: Parantheses are needed
        }
        fe[row] = detA/6.0/4.0*f_func(xm, ym, zm);
    }
}


void P1_3d::CalcLaplace(
    int const ial[4], double const xc[], 
    vector<vector<double>> &ske, vector<double> &fe) const
{
    CalcLaplace(ial,xc,ske,fe,rhs_lap3);
}



void P1_2vec_3d::CalcLaplace(
    int const ial[34], double const xc[], 
    vector<vector<double>> &ske, vector<double> &fe) const
{
    cout << "P1_2vec_3d::CalcLaplace\n";
    assert(nVertices_loc()==10);
    assert(nDOFs_loc()==34);
    //cout << ske.size() << "  " << nDOFs_loc() << endl;
    assert(static_cast<int>(ske.size())==nDOFs_loc());
    
    // Dummy filling for first test
    for (size_t row=0; row<ske.size(); ++row)
    {
        fill(ske[row].begin(),ske[row].end(),-1);
        ske[row][row] = ske.size()+0.5;
        fe[row] = 1.0/ske.size();
    }
    //   Here we have to fill ske, fe

}



