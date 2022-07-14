#pragma once

#include "geom.h"
#include "getmatrix.h"
#include "userset.h"
#include "utils.h"
#include "vdop.h"

#include <array>
#include <cassert>
#include <type_traits>
#include <variant>   // class variant
#include <vector>

//!< linear test functions for P1_3D element
static
std::array<std::function<double(double,double,double)>,4> const 
 phi_P1_3d{
          [](double xi_0,double xi_1,double xi_2) {return 1-xi_0-xi_1-xi_2;},
          [](double xi_0,double,double) {return xi_0;},
          [](double,double xi_1,double) {return xi_1;},
          [](double,double,double xi_2) {return xi_2;}
          };

//!< derivatives of linear test functions for P1_3D element
static
std::array<std::array<std::function<double(double,double,double)>,3>  ,4> const
  //dphi_P1_3d{
              //[](double,double,double) {return -1.0;}, [](double,double,double) {return -1.0;}, [](double,double,double) {return -1.0;} , 
              //[](double,double,double) {return  1.0;}, [](double,double,double) {return  0.0;}, [](double,double,double) {return  0.0;} , 
              //[](double,double,double) {return  0.0;}, [](double,double,double) {return  1.0;}, [](double,double,double) {return  0.0;} , 
              //[](double,double,double) {return  0.0;}, [](double,double,double) {return  0.0;}, [](double,double,double) {return  1.0;}  
          //};
// Always use double {{}} for initializer list for array, see https://stackoverflow.com/questions/22501368/why-wasnt-a-double-curly-braces-syntax-preferred-for-constructors-taking-a-std
  dphi_P1_3d{{
             {{ [](double,double,double) {return -1.0;}, [](double,double,double) {return -1.0;}, [](double,double,double) {return -1.0;} }}, 
             {{ [](double,double,double) {return  1.0;}, [](double,double,double) {return  0.0;}, [](double,double,double) {return  0.0;} }}, 
             {{ [](double,double,double) {return  0.0;}, [](double,double,double) {return  1.0;}, [](double,double,double) {return  0.0;} }}, 
             {{ [](double,double,double) {return  0.0;}, [](double,double,double) {return  0.0;}, [](double,double,double) {return  1.0;} }} 
          }};          

//!< quadratic test functions for P2_3D element
static
std::array<std::function<double(double,double,double)>,10> const 
 phi_P2_3d{
          [](double xi_0,double xi_1,double xi_2) {return phi_P1_3d[0](xi_0,xi_1,xi_2)*(2*phi_P1_3d[0](xi_0,xi_1,xi_2)-1);},
          [](double xi_0,double xi_1,double xi_2) {return phi_P1_3d[1](xi_0,xi_1,xi_2)*(2*phi_P1_3d[1](xi_0,xi_1,xi_2)-1);},
          [](double xi_0,double xi_1,double xi_2) {return phi_P1_3d[2](xi_0,xi_1,xi_2)*(2*phi_P1_3d[2](xi_0,xi_1,xi_2)-1);},
          [](double xi_0,double xi_1,double xi_2) {return phi_P1_3d[3](xi_0,xi_1,xi_2)*(2*phi_P1_3d[3](xi_0,xi_1,xi_2)-1);},
          [](double xi_0,double xi_1,double xi_2) {return 4*phi_P1_3d[0](xi_0,xi_1,xi_2)*phi_P1_3d[1](xi_0,xi_1,xi_2);},
          [](double xi_0,double xi_1,double xi_2) {return 4*phi_P1_3d[1](xi_0,xi_1,xi_2)*phi_P1_3d[2](xi_0,xi_1,xi_2);},
          [](double xi_0,double xi_1,double xi_2) {return 4*phi_P1_3d[0](xi_0,xi_1,xi_2)*phi_P1_3d[2](xi_0,xi_1,xi_2);},
          [](double xi_0,double xi_1,double xi_2) {return 4*phi_P1_3d[0](xi_0,xi_1,xi_2)*phi_P1_3d[3](xi_0,xi_1,xi_2);},
          [](double xi_0,double xi_1,double xi_2) {return 4*phi_P1_3d[1](xi_0,xi_1,xi_2)*phi_P1_3d[3](xi_0,xi_1,xi_2);},
          [](double xi_0,double xi_1,double xi_2) {return 4*phi_P1_3d[2](xi_0,xi_1,xi_2)*phi_P1_3d[3](xi_0,xi_1,xi_2);},
          };


// static polymorphism ?
// yes, via template argument for FEM_Matrix_2<ELEM>
// R.Grimm: https://www.modernescpp.com/index.php/dynamic-and-static-polymorphism
/**
 * Basis class for finite elements.
 */ 
class Element
{
    public:
    constexpr Element(int ndim, int nvertices, int ndof)
     : _ndim(ndim), _nvertices(nvertices), _ndof(ndof)  {}
    
    virtual ~Element();
    //~Element() = default;

    /**
     * @return Space dimension.
     */    
    constexpr int nDim_loc() const
    {return _ndim;}
    
    /**
     * @return number of geometric vertices defining the element
     */    
    constexpr int nVertices_loc() const
    {return _nvertices;}

    /**
     * @return number of degrees of freedem in element
     */    
    constexpr int nDOFs_loc() const
    {return _ndof;}
                
    /**
     * Derives the global indices for the DOFs in each element from the 
     * global indices of geometric vertices @p ia_geom.
     * 
     * @param[in] ia_geom global geometric vertex indices[nelem*_nvertices]
     * 
     * @return global indices for the DOFs [nelem*_ndof]
     */  
    virtual 
    std::vector<int> getGlobalDOFIndices(std::vector<int> const &ia_geom) const = 0;
    
    /**
     * Allocates the local finite element matrix.
     * 
     * @return matrix [_ndof][_ndof]
     */  
    std::vector<std::vector<double>> createElemMatrix() const
    {
        return std::vector<std::vector<double>>(_ndof,std::vector<double>(_ndof,0.0));
    }

    /**
     * Allocates the local finite element load vector.
     * 
     * @return vector [_ndof]
     */  
    std::vector<double> createElemVector() const
    {
        return std::vector<double>(_ndof,0.0);
    }

    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     */        
    void CalcLaplace(
         int const ial[3], 
         double const xc[], 
         std::vector<std::vector<double>> &ske, 
         std::vector<double> &fe) const;

    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @p fe will be calculated according to function @p f_func (x,y,z)
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     * @param[out] f_func function @p f_func (x,y,z)
     */      
    void CalcLaplace(
         int const ial[3], double const xc[], 
         std::vector<std::vector<double>> &ske, std::vector<double> &fe,
         const std::function<double(double, double, double)> &f_func
         ) const;
    
    private:
    int const _ndim;      //!< space dimension {2,3}
    int const _nvertices; //!< number of geometric vertices
    int const _ndof;      //!< number of degrees of freedom;
    // TODO: Add polynomial degree?
};


/**
 *  Linear triangular P1-element.
 *  The local vertex numbering follows Jung/Langer, Tabelle 4.3 (page 256).
 */ 
class P1_2d: public Element
{
    public:
    /**
     * @param[in] ndof_v degrees of freedom per vertex
     */ 
    constexpr     
    P1_2d(int ndof_v=1)
     : Element(2,3,3*ndof_v) { assert(0 == nDOFs_loc() % nVertices_loc() ); }
     
    virtual ~P1_2d() override {}
     
    /**
     * Derives the global indices for the DOFs in each element from the 
     * global indices of geometric vertices @p ia_geom.
     * 
     * @param[in] ia_geom global geometric vertex indices[nelem*_nvertices]
     * 
     * @return global indices for the DOFs [nelem*_ndof]
     */ 
    std::vector<int> getGlobalDOFIndices(std::vector<int> const &ia_geom) const override;

    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     */    
    void CalcLaplace(
        int const ial[3], 
        double const xc[], 
        std::vector<std::vector<double>> &ske, 
        std::vector<double> &fe) const;

    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @p fe will be calculated according to function @p f_func (x,y)
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     * @param[out] f_func function @p f_func (x,y)
     */      
    void CalcLaplace(
         int const ial[3], double const xc[], 
         std::vector<std::vector<double>> &ske, std::vector<double> &fe,
         const std::function<double(double, double)> &f_func
         ) const;
};

///**
 //*  Linear triangular P1-element.
 //*  The local vertex numbering follows Jung/Langer, Tabelle 4.3 (page 256).
 //*/ 
//class P1_2d_1dof: public P1_2d
//{
    //public:
    ///**
     //* one degree of freedom per vertex
     //*/ 
    ////constexpr P1_2d_1dof()
     ////: Element(2,3,3*1) { }
    //P1_2d_1dof()
     //: P1_2d(1) { }
    
    //virtual ~P1_2d_1dof() override {}
//};

/**
 *  Linear tetrahedral P1-element.
 *  The local vertex numbering follows Jung/Langer, Tabelle 4.3 (page 256).
 */ 
class P1_3d: public Element
{
    public:
    /**
     * @param[in] ndof_v degrees of freedom per vertex
     */
    constexpr
    P1_3d(int ndof_v=1)
     : Element(3,4,4*ndof_v) { assert(0 == nDOFs_loc() % nVertices_loc() ); }
     
    virtual ~P1_3d() override {}
     
    /**
     * Derives the global indices for the DOFs in each element from the 
     * global indices of geometric vertices @p ia_geom.
     * 
     * @param[in] ia_geom global geometric vertex indices[nelem*_nvertices]
     * 
     * @return global indices for the DOFs [nelem*_ndof]
     */ 
    std::vector<int> getGlobalDOFIndices(std::vector<int> const &ia_geom) const override;
    
    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     */    
    void CalcLaplace(int const ial[4], double const xc[], std::vector<std::vector<double>> &ske, std::vector<double> &fe) const;

    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @p fe will be calculated according to function @p f_func (x,y,z)
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     * @param[out] f_func function @p f_func (x,y,z)
     */      
    void CalcLaplace(
         int const ial[4], double const xc[], 
         std::vector<std::vector<double>> &ske, std::vector<double> &fe,
         const std::function<double(double, double, double)> &f_func    ) const;    
};

/**
 *  Quadratic tetrahedral P2-element.
 */ 
class P2_3d: public Element
{
    public:
    /**
     * @param[in] ndof_v degrees of freedom per vertex
     */
    constexpr
    P2_3d(int ndof_v=1)
     : Element(3,10,10*ndof_v) { assert(0 == nDOFs_loc() % nVertices_loc() ); }
     
    virtual ~P2_3d() override {}

    /**
     * Derives the global indices for the DOFs in each element from the 
     * global indices of geometric vertices @p ia_geom.
     * 
     * @param[in] ia_geom global geometric vertex indices[nelem*_nvertices]
     * 
     * @return global indices for the DOFs [nelem*_ndof]
     */ 
    std::vector<int> getGlobalDOFIndices(std::vector<int> const &ia_geom) const override;
};


/**
 *  Combined P1-P2 tetrahedral element with vector values in P2 vertices.
 * 
 *  Local numbering of dofs: 4 scalar values (linear test function) followed 
 *  by 10 vector values (quadratic test function) consisting of 3 components. 
 *
 * 
 */ 
class P1_2vec_3d: public Element
{
    public:
    /**
     * @param[in] ndof_v  max. degrees of freedom per vertex
     */
    constexpr
    P1_2vec_3d(int ndof_v=3)
     : Element(3,10,4+10*ndof_v) {  }
     
    virtual ~P1_2vec_3d() override {}

    /**
     * Derives the global indices for the DOFs in each element from the 
     * global indices of geometric vertices @p ia_geom.
     * 
     * @param[in] ia_geom global geometric vertex indices[nelem*_nvertices]
     * 
     * @return global indices for the DOFs [nelem*_ndof]
     */ 
    std::vector<int> getGlobalDOFIndices(std::vector<int> const &ia_geom) const override;

    /**
     * Calculates the element stiffness matrix @p ske and the element load vector @p fe
     * of one triangular element with linear shape functions.
     * @param[in]	ial	node indices of the element vertices [_nvertices]
     * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
     * @param[out] ske	element stiffness matrix [_ndof][_ndof]
     * @param[out] fe	element load vector [_ndof]
     */    
    void CalcLaplace(int const ial[34], double const xc[], std::vector<std::vector<double>> &ske, std::vector<double> &fe) const;    
    
};

/**
 * Copies only the first @p col1 columns from matrix @p v2 [n* @p col2]  with @p col1<= @p col2.
 *  
 * @param[in] v2   data[n*col2]
 * @param[in] col1 #columns in output
 * @param[in] col2 #columns in input
 * @return    data[n*col1]
 */ 
std::vector<int> shrinkColumns(std::vector<int> const &v2, int col2, int col1);

/**
 * Merges data @p v1 and data @p into a new vector/matrix.
 *  
 * @param[in] v1   data[n*col1]
 * @param[in] col1 #columns in @p v1
 * @param[in] v2   data[n*col2]
 * @param[in] col2 #columns in @p v2
 * @return    data[n*(col1+col2)]
 */
std::vector<int> mergeColumns(std::vector<int> const &v1, int col1, std::vector<int> const &v2, int col2);

/**
 * An index vector @p vin for scalar values is transferred into 
 * an index vector for @p ndof_v DOFs per scalar.
 * An optional index @p offset is applied (useful for merging two index vectors).
 *  
 * @param[in] vin    index[n]
 * @param[in] ndof_v #DOFs per original index
 * @param[in] offset index offset for output indices
 * @return    index[n*ndof_v]
 * 
 * @see mergeColumns
 */
std::vector<int> indexVectorBlowUp(std::vector<int> const &vin, int ndof_v=1, int offset=0);


//######################################################################

/**
 * FEM Matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */
template <class ELEM=P1_3d()>
class FEM_Matrix_2: public CRS_Matrix
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
       explicit FEM_Matrix_2(Mesh const & mesh, ELEM const &elem=P1_3d());
       //explicit FEM_Matrix_2(Mesh const & mesh);

       FEM_Matrix_2(FEM_Matrix_2 const &) = default;

      /**
        * Destructor.
        */
       ~FEM_Matrix_2() override;
       
    /**
     * Read DOF connectivity information (g1,g2,g3)_i.
     * @return connectivity vector [nelems*ndofs].
     */
    [[nodiscard]] const std::vector<int>  &GetConnectivity() const
    {
        return _ia_dof;
    }
    
    /**
     * Read DOF connectivity information (g1,g2,g3)_i.
     * @return connectivity vector [nelems*ndofs].
     */
    [[nodiscard]] const std::vector<int>  &GetConnectivityGeom() const
    {
        return _mesh.GetConnectivity();
    }
    

    /**
     * Determines all node to node connections from the dof based mesh.
      *
      * Faster than @p Node2NodeGraph_1().
      *
     * @return vector[k][] containing all connections of dof k, including to itself.
     */    
    std::vector<std::vector<int>> Node2NodeGraph() const
    {
        //std::cout << Nelems() << " " << NdofsElement() << " " << GetConnectivity().size() << std::endl;
        return ::Node2NodeGraph(Nelems(),NdofsElement(),GetConnectivity());
    }

    
    /**
     * Checks whether the mesh is compatible with the choosen finite element type.
     * @return true iff compatible.
     */
    bool CheckCompatibility() const;
    
    constexpr int Nelems() const
    { return _mesh.Nelems();}
    
    constexpr int Nnodes() const
    { return _mesh.Nnodes();}
    
    constexpr int NdofsElement() const
    { return _elem.nDOFs_loc(); }
    
    constexpr int NverticesElement() const
    { return _elem.nVertices_loc(); }
    
    constexpr int Ndims() const
    { return _elem.nDim_loc(); }    

       /**
        * Generates the sparse matrix pattern and overwrites the existing pattern.
        *
        * The sparse matrix pattern is generated but the values are 0.
       */
       void Derive_Matrix_Pattern();
       //{
           ////Derive_Matrix_Pattern_slow();
           //Derive_Matrix_Pattern_fast();
           //CheckRowSum();
       //}
       //void Derive_Matrix_Pattern_fast();
       //void Derive_Matrix_Pattern_slow();


        /**
        * Calculates the entries of f.e. stiffness matrix for the Laplace operator
        * and load/rhs vector @p f.
        * No memory is allocated.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        * @warning Only linear elements (P1_2d, P1_3d) are supported.
        * @see P1_2d
        * @see P1_3d
        */
       void CalculateLaplace(std::vector<double> &f);
       
        /**
        * Calculates the entries of f.e. stiffness matrix for the Laplace operator in 2D
        * and load/rhs vector @p f according to function @p f_func,
        *
        * @param[in,out] f (preallocated) rhs/load vector
        * @param[in]     f_func function f(x,y)
        */
       void CalculateLaplace(
           std::vector<double> &f, 
           std::function<double(double, double)> const &f_func);

        /**
        * Calculates the entries of f.e. stiffness matrix for the Laplace operator in 3D
        * and load/rhs vector @p f according to function @p f_func,
        *
        * @param[in,out] f (preallocated) rhs/load vector
        * @param[in]     f_func function f(x,y,z)
        */
       void CalculateLaplace(
           std::vector<double> &f, 
           std::function<double(double, double, double)> const &f_func);

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
       //void AddElem_3(int const ial[3], double const ske[3][3], double const fe[3], std::vector<double> &f);
       void AddElem(
           int const ial[], 
           std::vector<std::vector<double>> const& ske, 
           std::vector<double> const &fe, 
           std::vector<double> &f);

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
        * Applies Dirichlet boundary conditions to stiffness matrix and to load vector @p f.
        * The <a href="https://www.jstor.org/stable/2005611?seq=1#metadata_info_tab_contents">penalty method</a>
        * is used for incorporating the given values @p u 
        * at the surface of the bounding box [@p xl, @p xh]x[@p yl, @p yh]x[@p zl, @p zh].
        * 
        * Skipping @p zl and @p zh reduced the box to a rectangle in 2D.
        *
        * @param[in]     u  (global) vector with Dirichlet data
        * @param[in,out] f  load vector
        * @param[in]     xl lower  value x-bounds
        * @param[in]     xh higher value x-bounds
        * @param[in]     yl lower  value y-bounds
        * @param[in]     yh higher value y-bounds
        * @param[in]     zl lower  value z-bounds
        * @param[in]     zh higher value z-bounds
        */
       void ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            double xl, double xh, double yl, double yh, 
            double zl=-std::numeric_limits<double>::infinity(), 
            double zh= std::numeric_limits<double>::infinity() );

    private:
       Mesh const     & _mesh;   //!< discretization (contains element connectivity regarding geometric vertices)
       ELEM const       _elem;   //!< element type
       std::vector<int> _ia_dof; //!< element connectivity regarding DOFs in element

};

//######################################################################

template <class ELEM>
FEM_Matrix_2<ELEM>::FEM_Matrix_2(Mesh const &mesh, ELEM const &elem)
//FEM_Matrix_2<ELEM>::FEM_Matrix_2(Mesh const &mesh)
    : CRS_Matrix(), _mesh(mesh), _elem(elem), _ia_dof()
    //: CRS_Matrix(), _mesh(mesh), _elem(ELEM()), _ia_dof()
{
    assert(CheckCompatibility());
    
    _ia_dof = _elem.getGlobalDOFIndices(GetConnectivityGeom());
    Derive_Matrix_Pattern();                      // requires _ia_dof
    return;
}

template <class ELEM>
FEM_Matrix_2<ELEM>::~FEM_Matrix_2()
{}

template <class ELEM>
bool FEM_Matrix_2<ELEM>::CheckCompatibility() const
{
    int const ndim_fem  = _elem.nDim_loc();
    int const nDOF_fem  = _elem.nDOFs_loc();
    int const nvert_fem = _elem.nVertices_loc();
    
    int const ndim = _mesh.Ndims();
    int const nDOF = _mesh.NdofsElement();
    int const nvert= _mesh.NverticesElement();
    std::cout << "##  FEM_Matrix_2::CheckCompatibility()  ##" << std::endl;
    
    bool ok{true};
    if (ndim_fem!=ndim)
    {
        std::cout << "Space dimensions mismatch: " << ndim << " vs. " << ndim_fem << std::endl;
        ok = false;
    }
    if (ok && nvert_fem!=nvert)
    {
        std::cout << "#local vertices mismatch: " << nvert << " vs. " << nvert_fem << std::endl;
        ok = false;
    }
    if (ok && nDOF_fem < nDOF)
    {
        std::cout << "#local DOFs to small: " << nDOF << " vs. " << nDOF_fem << std::endl;
        ok = false;
    }
    
    return true;
}

template <class ELEM>
void FEM_Matrix_2<ELEM>::Derive_Matrix_Pattern()
{
    if (_ia_dof.empty())
    {
        _ia_dof = _elem.getGlobalDOFIndices(GetConnectivityGeom());
    }
    std::cout << "\n############   FEM_Matrix::Derive_Matrix_Pattern ";
    MyTimer tstart;    //tstart.tic();
    
    //int const nelem(Nelems());
    //int const ndof_e(NdofsElement());
    auto const &ia(GetConnectivity());
    
    //std::cout << ia << std::endl;
    
//  Determine the number of matrix rows
    //_nrows = *max_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem);
    _nrows = *max_element(ia.cbegin(), ia.cend()); // GH??
    ++_nrows;                                 // node numberng: 0 ... nnode-1
    //assert(*min_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem) == 0); // numbering starts with 0 ?
    assert(*min_element(ia.cbegin(), ia.cend()) == 0); // numbering starts with 0 ? // GH??

// CSR data allocation
    _id.resize(_nrows + 1);                  // Allocate memory for CSR row pointer
//##########################################################################
    //auto const v2v = _mesh.Node2NodeGraph();               // GH TODO:  Hier Node2NodeGraph mit _ia_dof
    auto const v2v = Node2NodeGraph();               // GH TODO:  Hier Node2NodeGraph mit _ia_dof
    _nnz = 0;                         // number of connections
    _id[0] = 0;                       // start of matrix row zero
    for (size_t v = 0; v < v2v.size(); ++v ) {
        _id[v + 1] = _id[v] + v2v[v].size();
        _nnz += v2v[v].size();
    }
    assert(_nnz == _id[_nrows]);             // GH: TODO
    _sk.resize(_nnz,-12345.0);               // Allocate memory for CSR values vector

// CSR data allocation
    _ik.resize(_nnz);                        // Allocate memory for CSR column index vector
// Copy column indices
    int kk = 0;
    for (const auto & v : v2v) {
        for (size_t vi = 0; vi < v.size(); ++vi) {
            _ik[kk] = v[vi];
            ++kk;
        }
    }
    _ncols = *max_element(_ik.cbegin(), _ik.cend());  // maximal column number
    ++_ncols;                                         // node numbering: 0 ... nnode-1
    //cout << _nrows << "  " << _ncols << endl;
    assert(_ncols == _nrows);

    std::cout << "finished in  " <<  tstart.toc()  << " sec.    ########\n";

    return;
}

template <class ELEM>
void FEM_Matrix_2<ELEM>::AddElem(int const ial[], std::vector<std::vector<double>> const &ske, std::vector<double> const &fe, std::vector<double> &f)
{
    //assert(NdofsElement() == 3);               // only for triangular, linear elements
    for (int i = 0; i < NdofsElement(); ++i) {
        const int ii  = ial[i];           // row ii (global index)
        for (int j = 0; j < NdofsElement(); ++j) {     // no symmetry assumed
            const int jj = ial[j];        // column jj (global index)
            const int ip = fetch(ii, jj);       // find column entry jj in row ii
#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
            if (ip < 0) {      // no entry found !!
                std::cout << "Error in AddElem: (" << ii << "," << jj << ") ["
                          << ial[0] << "," << ial[1] << "," << ial[2] << "]\n";
                assert(ip >= 0);
            }
#endif
            #pragma omp atomic
            _sk[ip] += ske[i][j];
        }
        #pragma omp atomic
        f[ii] += fe[i];
    }
}

template <class ELEM>
void FEM_Matrix_2<ELEM>::CalculateLaplace(
    std::vector<double> &f, 
    std::function<double(double, double)> const &f_func)          // 2D
{
    assert(_elem.nDim_loc()==2);                                  // 2D
    std::cout << "\n############   FEM_Matrix::CalculateLaplace f_func 2D";
    //std::cout << _nnz << " vs. " << _id[_nrows] << "  " << _nrows<< std::endl;
    assert(_nnz == _id[_nrows]);
    
    MyTimer tstart;    //tstart.tic();
    fill(_sk.begin(),_sk.end(),0.0);              // set     matrix entries to zero
    fill(  f.begin(),  f.end(),0.0);              // set rhs vector entries to zero

    auto const nelem    = Nelems();
    auto const nvert_e  = NverticesElement();     // geom. vertices per element
    auto const ndof_e   = NdofsElement();         // DOFs per element
    auto const &ia_geom = GetConnectivityGeom();  // geometric connectivity
    auto const &ia_crs  = GetConnectivity();      // DOF connectivity in CRS matrix
    auto const &xc      = _mesh.GetCoords();      // Coordinates

    #pragma omp parallel
    {
    auto ske(_elem.createElemMatrix());
    auto fe (_elem.createElemVector());
    //#pragma omp parallel for private(ske,fe)    // GH: doesn't work correctly with vector<vector<double>>
    #pragma omp for    
    for (int i = 0; i < nelem; ++i) {             //  Loop over all elements
        _elem.CalcLaplace(ia_geom.data() + nvert_e * i, xc.data(), ske, fe,f_func);
        AddElem(ia_crs.data() + ndof_e * i, ske, fe, f);
    }
    }

    std::cout << "finished in  " <<  tstart.toc()  << " sec.    ########\n";
    //Debug();
    return;
}

template <class ELEM>
void FEM_Matrix_2<ELEM>::CalculateLaplace(
    std::vector<double> &f, 
    std::function<double(double, double, double)> const &f_func)  // 3D
    //const std::variant<std::function<double(double, double)>, std::function<double(double, double, double)>> &f_func)
{
    assert(_elem.nDim_loc()==3);                                  // 3D
    std::cout << "\n############   FEM_Matrix::CalculateLaplace f_func 3D";
    //std::cout << _nnz << " vs. " << _id[_nrows] << "  " << _nrows<< std::endl;
    assert(_nnz == _id[_nrows]);
    
    MyTimer tstart;    //tstart.tic();
    fill(_sk.begin(),_sk.end(),0.0);              // set     matrix entries to zero
    fill(  f.begin(),  f.end(),0.0);              // set rhs vector entries to zero

    auto const nelem    = Nelems();
    auto const nvert_e  = NverticesElement();     // geom. vertices per element
    auto const ndof_e   = NdofsElement();         // DOFs per element
    auto const &ia_geom = GetConnectivityGeom();  // geometric connectivity
    auto const &ia_crs  = GetConnectivity();      // DOF connectivity in CRS matrix
    auto const &xc      = _mesh.GetCoords();      // Coordinates

    #pragma omp parallel
    {
    auto ske(_elem.createElemMatrix());
    auto fe (_elem.createElemVector());
    //#pragma omp parallel for private(ske,fe)    // GH: doesn't work correctly with vector<vector<double>>
    #pragma omp for    
    for (int i = 0; i < nelem; ++i) {             //  Loop over all elements
        _elem.CalcLaplace(ia_geom.data() + nvert_e * i, xc.data(), ske, fe,f_func);
        AddElem(ia_crs.data() + ndof_e * i, ske, fe, f);
    }
    }

    std::cout << "finished in  " <<  tstart.toc()  << " sec.    ########\n";
    //Debug();
    return;
}

// only for linear elements and scalar functions
template <class ELEM>
void FEM_Matrix_2<ELEM>::CalculateLaplace(std::vector<double> &f)
{
    std::cout << "\n############   FEM_Matrix::CalculateLaplace ";
    if constexpr(std::is_same<ELEM,P1_3d>::value)
    {
        CalculateLaplace(f, rhs_lap3 );                           // 3D standard rhs
    }
    else if constexpr(std::is_same<ELEM,P1_2d>::value)
    { 
        CalculateLaplace(f, rhs_lap2 );                           // 2D standard rhs
    }
    else
    {
        assert(false);
    }
}


// 2D-part: copy from getmatrix.cpp

template <class ELEM>
void FEM_Matrix_2<ELEM>::ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f)
{
    if (2==_elem.nDim_loc())
    {
        auto const idx = _mesh.Index_DirichletNodes();      // GH: not available in 3D
        int const nidx = idx.size();

        for (int i = 0; i < nidx; ++i) {
            int const row = idx[i];
            for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
                int const col = _ik[ij];
                if (col == row) {
                    _sk[ij] = 1.0;
                    f[row]  = u[row];
                 }
                else {
                    int const id1 = fetch(col, row); // Find entry (col,row)
                    assert(id1 >= 0);
                    f[col] -= _sk[id1] * u[row];
                    _sk[id1] = 0.0;
                    _sk[ij]  = 0.0;
                }
            }
        }
    }
    else
    {
        std::cout << "FEM_Matrix_2<ELEM>::ApplyDirichletBC not implemented in 3D" << std::endl;
        assert(false);
    }

    return;
}

template <class ELEM>
void FEM_Matrix_2<ELEM>::ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            double xl, double xh, double yl, double yh, double zl, double zh )
{
    std::vector<int> idx;
    if (3==_mesh.Ndims())
    {
        idx = _mesh.Index_DirichletNodes_Box(xl, xh, yl, yh, zl, zh);
        //std::cout << "####  3D idx: " << idx << std::endl;
    }
    else if (2==_mesh.Ndims())
    {
        idx = _mesh.Index_DirichletNodes_Box(xl, xh, yl, yh);
        //std::cout << "####  2D idx: " << idx << std::endl;
    }
    else
    {
        assert(false);
    }
    
    double const PENALTY = 1e6;
    int const nidx = idx.size();

    for (int row=0; row<nidx; ++row)
    {
        int const k = idx[row];
		int const id1 = fetch(k, k); // Find diagonal entry of row
		assert(id1 >= 0);
		_sk[id1] += PENALTY;		// matrix weighted scaling feasible
        f[k] += PENALTY * u[k];
    }

    return;
}


