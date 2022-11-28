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
/*
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
    for (int d=0; d<ndof_v; ++d)
    {
		 int const kv=d*vin.size();
         for (size_t k=0; k<vin.size(); ++k)
         {
             wout.at(k+kv) = vin.at(k)+offset+k;
         }
     }
    
    return wout;
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
    int const nP1_dofs = *max_element(cbegin(ia_p1),cend(ia_p1)) +1;  // offset for next numbering
    _nP1_nodes=nP1_dofs;
    _nP2_nodes=*max_element(cbegin(ia_geom_P2),cend(ia_geom_P2)) +1;
    
    // second: P2 part - vector
    auto const ia_p2(indexVectorBlowUp(ia_geom_P2,3,nP1_dofs));
    
    std::cout << "ia_p1: "  <<  ia_p1 << endl;
    std::cout << "ia_p2: "  <<  ia_p2 << endl;
    std::cout << "ia_geom_P1: "  <<  ia_geom_P1 << endl;
    std::cout << "ia_geom_P2: "  <<  ia_geom_P2 << endl;
    // merging of  ia_p1 and ia_p2
    int const col1=nv1*1;
    int const col2=nv2*3;
    vector<int> ia_dofs(mergeColumns(ia_p1, col1, ia_p2, col2));
    
    return ia_dofs;
}
*/

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

vector<int> mergeColumns(vector<int> const &v1, int col1, vector<int> const &v2, int col2, vector<int> const &v3, int col3, vector<int> const &v4, int col4)
{
    assert( v2.size()%col2 == 0 );
    assert( v1.size()%col1 == 0 );
    assert( v1.size()/col1 == v2.size()/col2 );
    int const nrows = static_cast<int>((v2.size()+v3.size()+v4.size())/(col2+col3+col4));
    int const ncols = col1+col2+col3+col4;
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
        for (int d=0; d<col3; ++d)
        {
            vm.at(k*ncols+d+col1+col2) = v3.at(k*col3+d);  // copy data v2
        }
        for (int d=0; d<col4; ++d)
        {
            vm.at(k*ncols+d+col1+col2+col3) = v4.at(k*col4+d);  // copy data v2
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
    int const nP1_dofs = *max_element(cbegin(ia_p1),cend(ia_p1)) +1;  // offset for next numbering
    _nP1_nodes=nP1_dofs;
    int const nP2_dofs=*max_element(cbegin(ia_geom_P2),cend(ia_geom_P2)) +1;
    _nP2_nodes=nP2_dofs;
    
    // second: P2 part - vector
    auto const ia_p2(indexVectorBlowUp(ia_geom_P2,1,nP1_dofs));
    auto const ia_p3(indexVectorBlowUp(ia_geom_P2,1,nP1_dofs+nP2_dofs));
    auto const ia_p4(indexVectorBlowUp(ia_geom_P2,1,nP1_dofs+2*nP2_dofs));
    
    //std::cout << "ia_p1: "  <<  ia_p1 << endl;
    //std::cout << "ia_p2: "  <<  ia_p2 << endl;
    //std::cout << "ia_geom_P1: "  <<  ia_geom_P1 << endl;
    //std::cout << "ia_geom_P2: "  <<  ia_geom_P2 << endl;
    // merging of  ia_p1 and ia_p2
    int const col1=nv1*1;
    int const col2=nv2*1;
    int const col3=nv2*1;
    int const col4=nv2*1;
    
    vector<int> ia_dofs(mergeColumns(ia_p1, col1, ia_p2, col2, ia_p3, col3, ia_p4, col4));
    
    return ia_dofs;
}




////######################################################################
/*
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



void P1_3d:: CalcElem_heat_equation_crank_nichelson(int const ial[4], double const xc[], vector<vector<double>> &ske, std::vector<double> &fe, const std::vector<double> &u_old, 
const double dt, const double t, const double c) const
{
    double tau=dt/c;

    const int  i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x1 = xc[i2 + 0] - xc[i1 + 0],  y1 = xc[i2 + 1] - xc[i1 + 1], z1 = xc[i2 + 2] - xc[i1 + 2],
                 x2 = xc[i3 + 0] - xc[i1 + 0],  y2 = xc[i3 + 1] - xc[i1 + 1], z2 = xc[i3 + 2] - xc[i1 + 2],
                 x3 = xc[i4 + 0] - xc[i1 + 0],  y3 = xc[i4 + 1] - xc[i1 + 1], z3 = xc[i4 + 2] - xc[i1 + 2],
                 x4 = xc[i3 + 0] - xc[i2 + 0],  y4 = xc[i3 + 1] - xc[i2 + 1], z4 = xc[i3 + 2] - xc[i2 + 2],
                 x5 = xc[i2 + 0] - xc[i4 + 0],  y5 = xc[i2 + 1] - xc[i4 + 1], z5 = xc[i2 + 2] - xc[i4 + 2],
                 x6 = xc[i4 + 0] - xc[i3 + 0],  y6 = xc[i4 + 1] - xc[i3 + 1], z6 = xc[i4 + 2] - xc[i3 + 2];
    const double jac = fabs(x3*(y1*z2-y2*z1)+y3*(x2*z1-x1*z2)+z3*(x1*y2-y1*x2));

    ske[0][0] = (1/(6*jac)) * ((x4*x4+y4*y4+z4*z4)*(x5*x5+y5*y5+z5*z5)-(x4*x5+y4*y5+z4*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][1] = (1/(6*jac)) * ((x2*x6+y2*y6+z2*z6)*(x4*x6+y4*y6+z4*z6)-(x2*x4+y2*y4+z2*z4)*(x6*x6+y6*y6+z6*z6));
    ske[0][2] = (1/(6*jac)) * ((x3*x4+y3*y4+z3*z4)*(x5*x5+y5*y5+z5*z5)-(x3*x5+y3*y5+z3*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][3] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x4*x6+y4*y6+z4*z6));
    ske[1][0] = ske[0][1];
    ske[1][1] = (1/(6*jac)) * ((x2*x2+y2*y2+z2*z2)*(x3*x3+y3*y3+z3*z3)-(x2*x3+y2*y3+z2*z3)*(x2*x3+y2*y3+z2*z3));
    ske[1][2] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x3*x6+y3*y6+z3*z6));
    ske[1][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x4*x6+y4*y6+z4*z6)-(x1*x6+y1*y6+z1*z6)*(x3*x4+y3*y4+z3*z4));
    ske[2][0] = ske[0][2];
    ske[2][1] = ske[1][2];
    ske[2][2] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x1*x3+y1*y3+z1*z3));
    ske[2][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x1*x4+y1*y4+z1*z4)-(x1*x1+y1*y1+z1*z1)*(x3*x4+y3*y4+z3*z4));
    ske[3][0] = ske[0][3];
    ske[3][1] = ske[1][3];
    ske[3][2] = ske[2][3];
    ske[3][3] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x1*x4+y1*y4+z1*z4));
    

    const double xm    = (xc[i1 + 0] + xc[i2 + 0] + xc[i3 + 0]+ xc[i4 + 0]) / 4.0,
                 ym    = (xc[i1 + 1] + xc[i2 + 1] + xc[i3 + 1]+ xc[i4 + 1]) / 4.0,
                 zm    = (xc[i1 + 2] + xc[i2 + 2] + xc[i3 + 2]+ xc[i4 + 2]) / 4.0;
    //fe[0] = fe[1] = fe[2] = 0.5 * jac * FunctF(xm, ym) / 3.0;
    fe[0] = fe[1] = fe[2]= fe[3] =  0*(jac/6.0) * fNice(t+dt/2, xm, ym, zm) / 4.0; // jac*0.5/3 is integral of one hat function over the triangle (volume of pyramid)

    // add contributions from crank nichelson to right hand-side

    fe[0] += -ske[0][0]*tau/2*u_old.at(ial[0])-ske[0][1]*tau/2*u_old.at(ial[1])-ske[0][2]*tau/2*u_old.at(ial[2])-ske[0][3]*tau/2*u_old.at(ial[3]);
    fe[1] += -ske[1][0]*tau/2*u_old.at(ial[0])-ske[1][1]*tau/2*u_old.at(ial[1])-ske[1][2]*tau/2*u_old.at(ial[2])-ske[1][3]*tau/2*u_old.at(ial[3]);
    fe[2] += -ske[2][0]*tau/2*u_old.at(ial[0])-ske[2][1]*tau/2*u_old.at(ial[1])-ske[2][2]*tau/2*u_old.at(ial[2])-ske[2][3]*tau/2*u_old.at(ial[3]);
    fe[3] += -ske[3][0]*tau/2*u_old.at(ial[0])-ske[3][1]*tau/2*u_old.at(ial[1])-ske[3][2]*tau/2*u_old.at(ial[2])-ske[3][3]*tau/2*u_old.at(ial[3]);


    ske[0][0] = ske[0][0]*tau/2 + jac/12.0;
    ske[0][1] = ske[0][1]*tau/2 + jac/24.0;
    ske[0][2] = ske[0][2]*tau/2 + jac/24.0;
    ske[0][3] = ske[0][2]*tau/2 + jac/24.0;
    ske[1][0] = ske[1][0]*tau/2 + jac/24.0;
    ske[1][1] = ske[1][1]*tau/2 + jac/12.0;
    ske[1][2] = ske[1][2]*tau/2 + jac/24.0;
    ske[1][3] = ske[1][2]*tau/2 + jac/24.0;
    ske[2][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[2][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[2][2] = ske[2][2]*tau/2 + jac/12.0;
    ske[2][3] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[3][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[3][2] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][3] = ske[2][2]*tau/2 + jac/12.0;


    // add contributions from mass matrix/time derivative to rhs

    fe[0] += jac/12.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[1] += jac/24.0*u_old.at(ial[0]) + jac/12.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[2] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/12.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[3] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/12.0*u_old.at(ial[3]);

}

*/






void P1_2vec_3d::CalcElem_Navier_Stokes(
    int const ial[10], double const xc[], std::vector<std::vector<double>> &ske, std::vector<double> &fe, const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
    const double dt, const double mu, const double lambda, const double kp, const double t_ni) const     
{
    //cout << "P1_2vec_3d::CalcLaplace\n";
    assert(nVertices_loc()==10);
    assert(nDOFs_loc()==34);
    //cout << ske.size() << "  " << nDOFs_loc() << endl;
    assert(static_cast<int>(ske.size())==nDOFs_loc());
    
    // Zero filling for first test
//    for (size_t row=0; row<ske.size(); ++row)
//    {
//        fill(ske[row].begin(),ske[row].end(),0.0);
//        fe[row] = 0.0;
//    }
    
    ::CalcElem_Navier_Stokes(ial, xc, ske, fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, mu, lambda, kp, t_ni);
      
}


/*
void P1_2vec_3d::CalcLaplace(
    int const ial[10], double const xc[], 
    vector<vector<double>> &ske, vector<double> &fe,
    const function<double(double, double, double)> &f_func ) const   
{
    cout << "P1_2vec_3d::CalcLaplace\n";
    assert(nVertices_loc()==10);
    assert(nDOFs_loc()==34);
    //cout << ske.size() << "  " << nDOFs_loc() << endl;
    assert(static_cast<int>(ske.size())==nDOFs_loc());
    
    // Zero filling for first test
    for (size_t row=0; row<ske.size(); ++row)
    {
        fill(ske[row].begin(),ske[row].end(),0.0);
        fe[row] = 0.0;
    }
    
        {
		double sa_ske[4][4];
		double sa_fe[4];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A00(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A00(..., sa_ske,sa_fe);
		int off_i = 0 + (1-1)*10;
        int off_j = 0 + (1-1)*10;
        for(int i=0; i<=3; ++i)
        {
			for(int j=0; j<=3; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			fe[off_i+i] = sa_fe[i];
		}
    }
    
    {
		double sa_ske[4][10];
		double sa_fe[4];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A01(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A01(..., sa_ske,sa_fe);
		int off_i = 0 + (1-1)*10;
        int off_j = 4 + (1-1)*10;
        for(int i=0; i<=3; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[4][10];
		double sa_fe[4];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A02(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A02(..., sa_ske,sa_fe);
		int off_i = 0 + (1-1)*10;
        int off_j = 4 + (2-1)*10;
        for(int i=0; i<=3; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[4][10];
		double sa_fe[4];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A03(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A03(..., sa_ske,sa_fe);
		int off_i = 0 + (1-1)*10;
        int off_j = 4 + (3-1)*10;
        for(int i=0; i<=3; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][4];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A10(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A10(..., sa_ske,sa_fe);
		int off_i = 4 + (1-1)*10;
        int off_j = 0 + (1-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=3; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			fe[off_i+i] = sa_fe[i];
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A11(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A11(..., sa_ske);
		int off_i = 4 + (1-1)*10;
        int off_j = 4 + (1-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A12(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A12(..., sa_ske);
		int off_i = 4 + (1-1)*10;
        int off_j = 4 + (2-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A13(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A13(..., sa_ske);
		int off_i = 4 + (1-1)*10;
        int off_j = 4 + (3-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][4];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A20(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A20(..., sa_ske,sa_fe);
		int off_i = 4 + (2-1)*10;
        int off_j = 0 + (1-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=3; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			fe[off_i+i] = sa_fe[i];
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A21(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A21(..., sa_ske);
		int off_i = 4 + (2-1)*10;
        int off_j = 4 + (1-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A22(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A22(..., sa_ske);
		int off_i = 4 + (2-1)*10;
        int off_j = 4 + (2-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A23(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A23(..., sa_ske);
		int off_i = 4 + (2-1)*10;
        int off_j = 4 + (3-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][4];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A30(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A30(..., sa_ske,sa_fe);
		int off_i = 4 + (3-1)*10;
        int off_j = 0 + (1-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=3; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			fe[off_i+i] = sa_fe[i];
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A31(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A31(..., sa_ske);
		int off_i = 4 + (3-1)*10;
        int off_j = 4 + (1-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
    
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A32(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A32(..., sa_ske);
		int off_i = 4 + (3-1)*10;
        int off_j = 4 + (2-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }
  
    {
		double sa_ske[10][10];
		double sa_fe[10];
		const std::vector<double> r_old_n;
		const std::vector<double> r_old_m;
		const std::vector<double> u_old_n;
		const std::vector<double> u_old_m;
		const std::vector<double> v_old_n;
		const std::vector<double> v_old_m;
		const std::vector<double> w_old_n;
		const std::vector<double> w_old_m;
		double dt, t, mu, lambda, kp;
		CalcElem_Navier_Stokes_A33(ial, xc, sa_ske, sa_fe, r_old_n, r_old_m, u_old_n, u_old_m, v_old_n, v_old_m, w_old_n, w_old_m, dt, t, mu, lambda, kp);
		//CalcElem_Navier_Stokes_A33(..., sa_ske);
		int off_i = 4 + (3-1)*10;
        int off_j = 4 + (3-1)*10;
        for(int i=0; i<=9; ++i)
        {
			for(int j=0; j<=9; ++j)
			{
				ske[off_i+i][off_j+j] = sa_ske[i][j];
			}
			//fe[off_i+i] = 0.0;
		}
    }


}
*/
