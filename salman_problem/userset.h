#ifndef USERSET_FILE
#define USERSET_FILE
#include <cmath>
/**
 * User function: f(@p x,@p y)
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value for right hand side f(@p x,@p y)
 */
double FunctF(double const x, double const y);

/**
 * User function: u(@p x,@p y)
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value for solution vector u(@p x,@p y)
 */
double FunctU(double const x, double const y);


/**
 * User function: f(@p x,@p y) = @f$ x^2 \sin(2.5\pi y)@f$.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value f(@p x,@p y)
 */
inline
double fNice(double const x, double const y)
{
    return x * x * std::sin(2.5 * 3.14159 * y);
}

/**
 * User function: f(@p x,@p y) = 0$.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value 0
 */
inline
double f_zero(double const x, double const y)
//double f_zero(double const  /*x*/, double const /*y*/)
{
    return 0.0 + 0.0*(x+y);
}

#endif
