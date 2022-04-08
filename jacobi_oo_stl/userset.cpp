#include "userset.h"
#include <cmath>


double FunctF(double const x , double const y)
{
// return  std::sin(3.14159*1*x)*std::sin(3.14159*1*y);
//  return 16.0*1024. ;
// return (double)1.0 ;
    return x * x * std::sin(2.5 * 3.14159 * y);
}

double FunctU(const double /* x */, double const /* y */)
{
    return 1.0 ;
}
