#ifndef JACSOLVE_FILE
#define JACSOLVE_FILE
#include "geom3.h"
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

void JacobiSolve(CRS_Matrix1 const &SK1, std::vector<double> const &f, std::vector<double> &u);


/**
 * Solves linear system of equations  K @p u = @p f  via the Jacobi iteration.
 * We use a distributed symmetric  CSR matrix @p SK and initial guess of the
 * solution is set to 0.
 *
 * In each smoothing step: @f$ \widehat{u} := u + \omega D^{-1}\left({f-K*u}\right) @f$
 *
 * @param[in]  SK	CSR matrix
 * @param[in]  f	distributed local vector storing the right hand side
 * @param[out] u	accumulated local vector storing the solution.
 * @param[in,out] r	auxiliary local vector.
 * @param[in] nsmooth	number of smoothing steps.
 * @param[in] omega	relaxation parameter.
 * @param[in] zero	initial solution @p u is assumed to be zero.
*/
void JacobiSmoother(Matrix const &SK1, std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> & r, int nsmooth=1, double const omega=1.0, bool zero=false);

/**
 * @brief Simple diagonale preconditioning.
 *
 * The residuum @p r scaled by the inverse diagonal of matríx @p SK results in the correction @p w.
 *
 * @f$ w :=  \omega D^{-1}*r @f$
 *
 * @param[in]  SK	matrix
 * @param[in]  r	distributed local vector storing the residuum
 * @param[out] w	accumulated local vector storing the correction.
 * @param[in] omega	relaxation parameter.
*/
void DiagPrecond(Matrix const &SK1, std::vector<double> const &r, std::vector<double> &w,
                 double const omega=1.0);



/**
 * @brief The Multigrid hierarchy including meshes, vectors and matrices, prolongations is stored.
*/




#endif
