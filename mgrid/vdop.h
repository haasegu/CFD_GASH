#pragma once
#include "utils.h"

#include <iostream>
#include <tuple>
#include <vector>

/** @brief  Element-wise vector divison x_k = y_k/z_k.
 *
 * @param[out] x  target vector
 * @param[in]  y  source vector
 * @param[in]  z  source vector
 *
 */
void vddiv(std::vector<double> & x, std::vector<double> const& y,
                                    std::vector<double> const& z);

/** @brief  Element-wise daxpy operation x(k) = y(k) + alpha*z(k).
 *
 * @param[out] x  target vector
 * @param[in]  y  source vector
 * @param[in]  alpha  scalar
 * @param[in]  z  source vector
 *
 */
void vdaxpy(std::vector<double> & x, std::vector<double> const& y,
                       double alpha, std::vector<double> const& z );


/** @brief  Calculates the Euclidean inner product of two vectors.
 *
 * @param[in]  x vector
 * @param[in]  y vector
 * @return Euclidean inner product @f$\langle x,y \rangle@f$
 *
 */
double dscapr(std::vector<double> const& x, std::vector<double> const& y);


inline
double L2_scapr(std::vector<double> const& x, std::vector<double> const& y)
{
    return dscapr(x,y)/static_cast<double>(x.size());
}

/**
 * Print entries of a vector.
 * @param[in] v	    vector values
*/
void DebugVector(std::vector<double> const &v);

/** @brief  Compares an STL vector with POD vector.
 *
 * The accuracy criteria @f$ |x_k-y_k| < \varepsilon \left({1+0.5(|x_k|+|y_k|)}\right) @f$
 * follows the book by
 * <a href="https://www.springer.com/la/book/9783319446592">Stoyan/Baran</a>, p.8.
 *
 * @param[in]  x    STL vector
 * @param[in]  n    length of POD vector
 * @param[in]  y    POD vector
 * @param[in]  eps  relative accuracy criteria (default := 0.0).
 * @return true iff pairwise vector elements are relatively close to each other.
 *
 */
bool CompareVectors(std::vector<double> const& x, int n, double const y[], double const eps=0.0);


/** Output operator for vector
 * 	@param[in,out] s	output stream, e.g. @p cout
 *  @param[in]     v    vector
 *
 *	@return    output stream
*/
template <class T>
std::ostream& operator<<(std::ostream &s, std::vector<T> const &v)
{
    for (auto vp: v)
    {
        s << vp << "  ";
    }
    return s;
}

/** Calculate the absolute difference between @p x and @p y.
 * 
 * 	@param[in]     x    vector
 *  @param[in]     y    vector
 *
 *	@return    absolute error vector
*/
std::vector<double> getAbsError(std::vector<double> const& x, std::vector<double> const& y);


/** Find the componentwise largest absolute difference between @p x and @p y.
 *  Optional via @p eps and @p nlarge the index and the value of the @p nlarge 
 *  error components will be printed to standard output.
 * 
 * 	@param[in]     x    vector
 *  @param[in]     y    vector
 *  @param[in]     eps  threshold
 *  @param[in]     nlarge  number of largest components to print.
 *
 *	@return    (max. error value, its index)
*/
std::tuple<double, int>  findLargestAbsError(
    std::vector<double> const& x, 
    std::vector<double> const& y,
    double eps=1e-6, int nlarge=0)           ;
    
    
bool vectorsAreEqual(
    std::vector<double> const& x, 
    std::vector<double> const& y,
    double eps=1e-6, int nlarge=0)           ;

