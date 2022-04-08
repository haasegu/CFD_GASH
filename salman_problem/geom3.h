#pragma once

#include "geom.h"

#include <iostream>
using namespace std;

/**
 * 2D finite element mesh of the square consiting of linear triangular elements.
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
	}
    
};
