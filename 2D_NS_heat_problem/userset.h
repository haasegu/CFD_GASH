#pragma once
#include <cmath>
/**
 * User function: f(@p x,@p y)
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value for right hand side f(@p x,@p y)
 */
double FunctF(double x, double y);

/**
 * User function: u(@p x,@p y)
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value for solution vector u(@p x,@p y)
 */
double FunctU(double x, double y);


/**
 * User function: f(@p x,@p y) = @f$ x^2 \sin(2.5\pi y)@f$.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value f(@p x,@p y)
 */
inline
double fNice(double const x, double const y)
{
    //return x * x * std::sin(2.5 * M_PI * y);               // solution u
    return std::sin(M_PI*2.5*y)*(M_PI*M_PI*2.5*2.5*x*x - 2); // -Laplacian(u)
}

/**
 * User function: f(@p x,@p y,@p z).
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @param[in] z		z-coordinate of discretization point
 * @return  value f(@p x,@p y,@p z)
 */
inline
double fNice(double const x, double const y, double const z)
{
    return std::sin(M_PI*2.5*y)*(M_PI*M_PI*2.5*2.5*x*x - 2)+z;
}

/**
 * User function: f(@p x,@p y) = @f$ x^2 \sin(2.5\pi y)@f$.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value f(@p x,@p y)
 */
inline
double fNice(const double t, double const x, double const y, double const z)
{
    return std::sin(M_PI*2.5*y*z)*(M_PI*M_PI*2.5*2.5*x*x - 2)*exp(-pow(t-1,2));
}

//Right hand side functions
inline
double fNice_x(double const x, double const y, double const z)
{
    return y*z+3*x*x*y*y*y*z-2-6*x*y;
}
inline
double fNice_y(double const x, double const y, double const z)
{
    return x*z+3*x*x*x*y*y*z-2*y*y;
}
inline
double fNice_z(double const x, double const y, double const z)
{
    return x*y+6*y*z+x*x*x*y*y*y;
}


inline
double u_boundary(double const x, double const y, double const z)
{
    //return x * x * std::sin(2.5 * M_PI * y);               // solution u
    return x*y*z;
    //sin(x*y*2*M_PI); // -Laplacian(u)
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

/**
 * User function: f(@p x,@p y) = 0$.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value 0
 */
inline
double f_zero(double const x, double const y, double const z)
//double f_zero(double const  /*x*/, double const /*y*/)
{
    return 0.0 + 0.0*(x+y+z);
}


// ---------------------------------------------------------------------

/**
 * solution u(@p x,@p y) of -Laplacian.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value u(@p x,@p y)
 * @see rhs_lap2
 */
inline
double sol_lap2(double const x, double const y)
{
    return x * x * std::sin(2.5 * M_PI * y);
}

/**
 * Right hand side f(@p x,@p y) of -Laplacian(u) with the appropriate function u.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @return  value f(@p x,@p y)
 * @see sol_lap2
 */
inline
double rhs_lap2(double const x, double const y)
{
    return std::sin(M_PI*2.5*y)*(M_PI*M_PI*2.5*2.5*x*x - 2);
}

// ---------------------------------------------------------------------

/**
 * solution u(@p x,@p y,@p z) of -Laplacian.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @param[in] z		z-coordinate of discretization point
 * @return  value u(@p x,@p y,@p z)
 * @see rhs_lap3
 */
inline
double sol_lap3(double const x, double const y, double const z)
{
    //return x*x*std::cos(M_PI*z)*std::sin(M_PI*y);
    return (x*x*(2*x - 3))*std::cos(M_PI*z)*std::cos(M_PI*y);
}

/**
 * Right hand side f(@p x,@p y,@p z) of -Laplacian(u) with the appropriate function u.
 * @param[in] x		x-coordinate of discretization point
 * @param[in] y		y-coordinate of discretization point
 * @param[in] z		z-coordinate of discretization point
 * @return  value f(@p x,@p y,@p z)
 * @see sol_lap3
 */
inline
double rhs_lap3(double const x, double const y, double const z)
{
    //return 2.0*std::cos(M_PI*z)*std::sin(M_PI*y)*(M_PI*M_PI*x*x - 1);
    return -2.0*std::cos(M_PI*z)*std::cos(M_PI*y)*(- 2*M_PI*M_PI*x*x*x + 3*M_PI*M_PI*x*x + 6*x - 3);
}

// ---------------------------------------------------------------------
