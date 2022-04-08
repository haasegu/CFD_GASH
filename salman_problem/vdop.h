#ifndef VDOP_FILE
#define VDOP_FILE
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


/** @brief  Calculates the Euclidian inner product of two vectors.
 *
 * @param[in]  x vector
 * @param[in]  y vector
 * @return Euclidian inner product @f$\langle x,y \rangle@f$
 *
 */
double dscapr(std::vector<double> const& x, std::vector<double> const& y);


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

#endif
