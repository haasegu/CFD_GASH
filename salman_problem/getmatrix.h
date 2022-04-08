#ifndef GETMATRIX_FILE
#define GETMATRIX_FILE

#include "geom.h"
#include <vector>

/**
 * Calculates the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions.
 * @param[in]	ial	node indices of the three element vertices
 * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coodinates of node k
 * @param[out] ske	element stiffness matrix
 * @param[out] fe	element load vector
 */
void CalcElem(int const ial[3], double const xc[], double ske[3][3], double fe[3]);

/**
 * Adds the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions to the appropriate positions in
 * the symmetric stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik)
 *
 * @param[in] ial   node indices of the three element vertices
 * @param[in] ske	element stiffness matrix
 * @param[in] fe	element load vector
 * @param[out] sk	vector non-zero entries of CSR matrix
 * @param[in] id	index vector containing the first entry in a CSR row
 * @param[in] ik	column index vector of CSR matrix
 * @param[out] f	distributed local vector storing the right hand side
 *
 * @warning Algorithm requires indices in connectivity @p ial in ascending order.
 *          Currently deprecated.
*/
void AddElem(int const ial[3], double const ske[3][3], double const fe[3],
             int const id[], int const ik[], double sk[], double f[]);


// #####################################################################
/**
 * Square matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */
class CRS_Matrix
{
    public:
       /**
        * Intializes the CRS matrix structure from the given discetization in @p mesh.
        *
        * The sparse matrix pattern is generated but the values are 0.
        *
        * @param[in] mesh given discretization
        *
        * @warning A reference to the discretization @p mesh is stored inside this class.
        *          Therefore, changing @p mesh outside requires also
        *          to call method @p Derive_Matrix_Pattern explicitely.
        *
        * @see Derive_Matrix_Pattern
        */
       explicit CRS_Matrix(Mesh const & mesh);

      /**
        * Destructor.
        */
       ~CRS_Matrix()
       {}

       /**
        * Generates the sparse matrix pattern and overwrites the existing pattern.
        *
        * The sparse matrix pattern is generated but the values are 0.
       */
       void Derive_Matrix_Pattern();

        /**
        * Calculates the entries of f.e. stiffness matrix and load/rhs vector @p f for the Laplace operator in 2D.
        * No memory is allocated.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       void CalculateLaplace(std::vector<double> &f);

       /**
        * Applies Dirichlet boundary conditions to stiffness matrix and to load vector @p f.
        * The <a href="https://www.jstor.org/stable/2005611?seq=1#metadata_info_tab_contents">penalty method</a>
        * is used for incorporating the given values @p u.
        *
        * @param[in]     u (global) vector with Dirichlet data
        * @param[in,out] f load vector
        */
       void ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f);

       /**
        * Extracts the diagonal elemenst of the sparse matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       void GetDiag(std::vector<double> &d) const;

       /**
        * Performs the matrix-vector product  w := K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     u vector
        */
       void Mult(std::vector<double> &w, std::vector<double> const &u) const;

        /**
        * Calculates the defect/residuum w := f - K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u vector
        */
       void Defect(std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const;

       /**
		 * Number rows in matrix.
		 * @return number of rows.
		 */
       int Nrows() const
          {return _nrows;}

       /**
		 * Show the matrix entries.
		 */
       void Debug() const;

       /**
		 * Finds in a CRS matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       int fetch(int row, int col) const;

      /**
        * Adds the element stiffness matrix @p ske and the element load vector @p fe
        * of one triangular element with linear shape functions to the appropriate positions in
        * the stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik).
        *
        * @param[in]     ial   node indices of the three element vertices
        * @param[in]     ske   element stiffness matrix
        * @param[in]     fe    element load vector
        * @param[in,out] f	   distributed local vector storing the right hand side
        *
        * @warning Algorithm assumes  linear triangular elements (ndof_e==3).
       */
       void AddElem_3(int const ial[3], double const ske[3][3], double const fe[3], std::vector<double> &f);

        /**
        * Compare @p this CRS matrix with an external CRS matrix stored in C-Style.
        *
        * The method prints statements on differences found.
        *
        * @param[in]     nnode  row number of external matrix
        * @param[in]     id     start indices of matrix rows of external matrix
        * @param[in]     ik     column indices of external matrix
        * @param[in]     sk     non-zero values of external matrix
        *
        * @return true iff all data are identical.
        */
       bool Compare2Old(int nnode, int const id[], int const ik[], double const sk[]) const;

    private:
       Mesh const & _mesh;      //!< reference to discretization
       int _nrows;              //!< number of rows in matrix
       int _nnz;                //!< number of non-zero entries
       std::vector<int> _id;    //!< start indices of matrix rows
       std::vector<int> _ik;    //!< column indices
       std::vector<double> _sk; //!< non-zero values

};


#endif
