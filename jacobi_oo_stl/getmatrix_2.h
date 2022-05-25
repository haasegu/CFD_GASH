#pragma once
#include "getmatrix.h"

/**
 * FEM Matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 * 
 * Based on quadratic tetrahedral elements.
 */
class FEM_Matrix_P2: public FEM_Matrix
{
    public:
       /**
        * Initializes the CRS matrix structure from the given discretization in @p mesh.
        *
        * The sparse matrix pattern is generated but the values are 0.
        *
        * @param[in] mesh given discretization
        *
        * @warning A reference to the discretization @p mesh is stored inside this class.
        *          Therefore, changing @p mesh outside requires also
        *          to call method @p Derive_Matrix_Pattern explicitly.
        *
        * @see Derive_Matrix_Pattern
        */
      
      explicit FEM_Matrix_P2(Mesh const & mesh, int ndof_v=10);
      
      
       FEM_Matrix_P2(FEM_Matrix_P2 const &) = default;

      /**
        * Destructor.
        */
       ~FEM_Matrix_P2() override;
       
       


        /**
        * Calculates the entries of f.e. stiffness matrix and load/rhs vector @p f for the Laplace operator in 2D.
        * No memory is allocated.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       //void CalculateLaplace(std::vector<double> &f);
       
       ///**
        //* Calculates the entries of f.e. matrix (stiffnes matrix + mass matrix using Crank-Nichelson implicit scheme) load/rhs vector @p f.
        //*
        //* @param[in,out] f (preallocated) rhs/load vector
        //*/
       //void CalculateLaplace_heat_equation(std::vector<double> &f, const std::vector<double> &u_old, const double dt, const double t, const double c);

       /**
        * Applies Dirichlet boundary conditions to stiffness matrix and to load vector @p f.
        * The <a href="https://www.jstor.org/stable/2005611?seq=1#metadata_info_tab_contents">penalty method</a>
        * is used for incorporating the given values @p u.
        *
        * @param[in]     u (global) vector with Dirichlet data
        * @param[in,out] f load vector
        */
       //void ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f);
       
       //void ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            //double xl, double xh, double yl, double yh,double zl, double zh );

       /**
        * Extracts the diagonal elements of the sparse matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
       */
       //void GetDiag(std::vector<double> &d) const;   // override in MPI parallel
       void GetDiag(std::vector<double> &d) const override   { std::cout << "GDDDDDD\n"; GetDiag_M(d); }
       
       // Solves non-M matrix problems for Jacobi iteration but not for MG
       //void GetDiag_M(std::vector<double> &d) const;
       
       //void Skalar2VectorMatrix(int ndof_v);
       

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
       //void AddElem_3(int const ial[4], double const ske[4][4], double const fe[4], std::vector<double> &f);
       
    /**
     * Global number of degrees of freedom (dof) for each finite element.
     * @return degrees of freedom per element.
     */
    //int NdofsVertex() const
    //{
        //return _ndof_v;
    //}

    private:
       //Mesh const & _mesh;      //!< reference to discretization
       //int  const _ndof_v;      //!< degrees of freedom per vertex (vector valued problems)

};
