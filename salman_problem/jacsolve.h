#ifndef JACSOLVE_FILE
#define JACSOLVE_FILE
#include "getmatrix.h"
#include <vector>


/**
 * Solves linear system of equations  K @p u = @p f  via the Jacobi iteration.
 * We use a distributed symmetric  CSR matrix @p SK and initial guess of the
 * solution is set to 0.
 * @param[in] SK	CSR matrix
 * @param[in] f		distributed local vector storing the right hand side
 * @param[out] u	accumulated local vector storing the solution.
*/
void JacobiSolve(CRS_Matrix const &SK, std::vector<double> const &f, std::vector<double> &u);


#endif
