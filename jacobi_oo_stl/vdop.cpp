#include "vdop.h"
#include <algorithm>
#include <cassert>               // assert()
#include <cmath>
#include <iostream>
#include <tuple>               // tuple
#include <vector>
using namespace std;


//void vddiv(size_t n, double x[], double const y[],
                              //double const z[])
//{
//#pragma omp parallel for default(none) shared(x,y,z) firstprivate(n)
    //for (size_t k = 0; k < n; ++k)
    //{
        //x[k] = y[k] / z[k];
    //}
    //return;
//}


void vddiv(vector<double> & x, vector<double> const& y,
                               vector<double> const& z)
{
    assert( x.size()==y.size() && y.size()==z.size() );
    size_t n = x.size();

#pragma omp parallel for    
    for (size_t k = 0; k < n; ++k)
    {
        x[k] = y[k] / z[k];
    }
    return;
}

//******************************************************************************

void vdaxpy(std::vector<double> & x, std::vector<double> const& y,
                       double alpha, std::vector<double> const& z )
{
    assert( x.size()==y.size() && y.size()==z.size() );
    size_t n = x.size();

#pragma omp parallel for    
    for (size_t k = 0; k < n; ++k)
    {
        x[k] = y[k] + alpha * z[k];
    }
    return;
}
//******************************************************************************

double dscapr(std::vector<double> const& x, std::vector<double> const& y)
{
    assert( x.size()==y.size());

    //return dscapr(x.size(),x.data(),y.data());

    size_t n = x.size();

    double    s = 0.0;
#pragma omp parallel for reduction(+:s)
    for (size_t k = 0; k < n; ++k)
    {
        s += x[k] * y[k];
    }

    return s;
}

//******************************************************************************
void DebugVector(vector<double> const &v)
{
    cout << "\nVector  (nnode = " << v.size() << ")\n";
    for (size_t j = 0; j < v.size(); ++j)
    {
        cout.setf(ios::right, ios::adjustfield);
        cout << v[j] << "   ";
    }
    cout << endl;

    return;
}
//******************************************************************************
bool CompareVectors(vector<double> const& x, int const n, double const y[], double const eps)
{
    bool bn = (static_cast<int>(x.size())==n);
    if (!bn)
    {
        cout << "#########   Error: " << "number of elements" << endl;
    }
    //bool bv = equal(x.cbegin(),x.cend(),y);
    bool bv = equal(x.cbegin(),x.cend(),y,
                          [eps](double a, double b) -> bool
                          { return std::abs(a-b)<eps*(1.0+0.5*(std::abs(a)+ std::abs(a))); }
    );
    if (!bv)
    {
        assert(static_cast<int>(x.size())==n);
        cout << "#########   Error: " << "values" << endl;
    }
    return bn && bv;
}

//******************************************************************************
vector<double> getAbsError(vector<double> const& x, vector<double> const& y)
{
    assert(size(x)==size(y));
    vector<double> err(size(x));
    
    for (size_t k=0; k<size(err); ++k)
    {
        //err[k] = std::abs( x[k]-y[k] );
        err[k] = x[k]-y[k];
    }
    return err;
}

//******************************************************************************
tuple<double, int>  findLargestAbsError(vector<double> const& x, vector<double> const& y,
    double eps, int nlarge)
{
    vector<double> err = getAbsError(x,y); 
    
    if (nlarge>0)
    {
        auto idx = sort_indexes_desc(err);
        for (int k=0; k<nlarge; ++k)
        {
            if ( err[idx[k]]>=eps )
            {
                cout << "err[" << idx[k] << "] = " << err[idx[k]] << endl;
            }
        }
    }
    
    auto ip = max_element(cbegin(err),cend(err));
    return make_tuple( *ip, ip-cbegin(err) );
}

bool vectorsAreEqual(vector<double> const& x, vector<double> const& y, double eps, int nlarge)
{
    bool bb=size(x)==size(y);
    bb += equal( cbegin(x),cend(x), cbegin(y),
    [eps](auto const a, auto const b)->bool 
    { return std::abs(a-b) <= eps*(1.0+(std::abs(a)+std::abs(b))/2.0); } 
    );
    
    if (!bb)
    {
        findLargestAbsError(x,y,eps,nlarge);
    }
    return bb;
}






//******************************************************************************

















