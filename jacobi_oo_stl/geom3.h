#pragma once

#include "geom.h"

#include <array>
#include <iostream>
#include <vector>
using namespace std;

/**
 * 3D finite element mesh linear tetrahedral elements.
 */
class Mesh_3d_4_matlab: public Mesh
{
public:
    /**
     * Reads mesh data from a binary file.
     *
     * File format, see ascii_write_mesh.m
     *
     * @param[in] fname file name
    */
    explicit Mesh_3d_4_matlab(std::string const &fname);
    
    
    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions.
     * @return index vector.
      *
      * @warning All boundary nodes are considered as Dirchlet nodes.
     */
    std::vector<int> Index_DirichletNodes() const override
    {
		cout <<   "Mesh_3d_4_matlab::Index_DirichletNodes()  not implemented" << endl;
        return std::vector<int>();        
	}
	
	std::vector<int> Index_DirichletNodes_Box(double xl, double xh, double yl, double yh,double zl, double zh );
};

/**
 * 3D finite element mesh linear tetrahedral elements (P1).
 * We start with a simple copy of Mesh_3d_4_matlab
 */
class Mesh_3d_P1_matlab: public Mesh
{
public:
    /**
     * Reads mesh data from a binary file.
     *
     * File format, see ascii_write_mesh.m
     *
     * @param[in] fname file name
    */
    explicit Mesh_3d_P1_matlab(std::string const &fname);
    
    
    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions.
     * @return index vector.
      *
      * @warning All boundary nodes are considered as Dirchlet nodes.
     */
    std::vector<int> Index_DirichletNodes() const override
    {
		cout <<   "Mesh_3d_P1_matlab::Index_DirichletNodes()  not implemented" << endl;
        return std::vector<int>();
	}
	
	std::vector<int> Index_DirichletNodes_Box(double xl, double xh, double yl, double yh,double zl, double zh );
};

/**
 * 3D finite element mesh quadratic tetrahedral elements (P2).
 */
class Mesh_3d_P2_matlab: public Mesh
{
public:
    /**
     * Reads mesh data from a binary file.
     *
     * File format, see ascii_write_mesh.m
     *
     * @param[in] fname file name
    */
    explicit Mesh_3d_P2_matlab(std::string const &fname);
    
   /**
     * Derives the local P2 mesh from the P1 mesh.
     * 
     * @param[in] p1 linear tetrahedral mesh (P1-element)
     * 
     */
    void deriveMeshFromP1(Mesh_3d_P1_matlab const &p1);

    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions.
     * @return index vector.
      *
      * @warning All boundary nodes are considered as Dirchlet nodes.
     */
    std::vector<int> Index_DirichletNodes() const override
    {
		cout <<   "Mesh_3d_P2_matlab::Index_DirichletNodes()  not implemented" << endl;
        return std::vector<int>();	
    }
	
	std::vector<int> Index_DirichletNodes_Box(double xl, double xh, double yl, double yh,double zl, double zh );
};

/**
 * Returns the vertex index of the arithmetic mean of vertices @p v1 and @p v2.
 * 
 * If that vertex is not already contained in the coordinate vector @xc then 
 * this new vertex is appended to @p xc.
 * 
 * @param[in]     v1    index of vertex
 * @param[in]     v2    index of vertex
 * @param[in,out] xc    coordinate vector [nnodes*ndim]
 * @param[in]     ndim  space dimension
 * @return vertex index of midpoint of vertices @p v1 and @p v2.
 * 
 */
int appendMidpoint(int v1, int v2, std::vector<double> &xc, int ndim=3);



/**
 * Determines the index of a vertex @p xm in the coordinate vector @p xc.
 * 
 * @param[in]     xm    one vertex
 * @param[in]     xc    vector of vertices [nnodes*ndim]
 * @param[in]     ndim  space dimension
 * @return index in vector or -1 in case the vertex is not contained in the vector.
 * 
 */
int getVertexIndex(std::vector<double> const &xm, std::vector<double> const &xc, int ndim=3);


/**
 * Compares two floating point numbers with respect to a sloppy accuracy.
 * 
 * @param[in]     a  number
 * @param[in]     b  number
 * @param[in]   eps  accuracy
 * @return @f$ |a-b| < \varepsilon @f$
 * 
 */
inline
bool equal(double a, double b, double eps=1e-6)
{
    return std::abs(b-a)<eps;
}









