#include "vdop.h"
#include <cassert>               // assert()
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;


void vddiv(vector<double> & x, vector<double> const& y,
                               vector<double> const& z)
{
    assert( x.size()==y.size() && y.size()==z.size() );
    size_t n = x.size();

#pragma omp parallel for default(none) shared(x,y,z) firstprivate(n)
    for (size_t k = 0; k < n; ++k)
    {
        x[k] = y[k] / z[k];
    }
    return;
}

//******************************************************************************

void vdaxpy(vector<double> & x, vector<double> const& y,
                  double alpha, vector<double> const& z )
{
    assert( x.size()==y.size() && y.size()==z.size() );
    size_t n = x.size();

#pragma omp parallel for default(none) shared(alpha, x,y,z)  firstprivate(n)
    for (size_t k = 0; k < n; ++k)
    {
        x[k] = y[k] + alpha * z[k];
    }
    return;
}
//******************************************************************************

double dscapr(vector<double> const& x, vector<double> const& y)
{
    assert( x.size()==y.size());
    size_t n = x.size();

    double    s = 0.0;
#pragma omp parallel for default(none) shared(x,y)  firstprivate(n) reduction(+:s)
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
    cout << endl;;

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
