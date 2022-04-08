
#include "geom3.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
using namespace std;

// ####################################################################

Mesh_3d_4_matlab::Mesh_3d_4_matlab(string const &fname)
 : Mesh(3,4,4)    // two dimensions, 3 vertices, 3 dofs
   //bedges(0)
{
    ifstream ifs(fname);
    if(!(ifs.is_open() && ifs.good()))
    {
        cerr << "Mesh_3d_4_matlab: Error cannot open file " << fname << endl;
        assert(ifs.is_open());
    }

    int const OFFSET(1);             // Matlab to C indexing
    cout << "ASCI file  " << fname << "  opened" << endl;

    // Read some mesh constants
    int nnode, ndim, nelem, nvert_e;
    ifs >> nnode >> ndim >> nelem >> nvert_e;
    cout << nnode << "  " << ndim << "  " << nelem << "  " << nvert_e <<endl;
    assert(ndim==3 && nvert_e==4);

    // Allocate memory
    Resize_Coords(nnode, ndim);                 // coordinates in 2D [nnode][ndim]
    Resize_Connectivity(nelem, nvert_e);        // connectivity matrix [nelem][nvert]

    // Read ccordinates
    auto &xc = GetCoords();
    for (int k = 0; k<nnode*ndim; ++k)
    {
        ifs >> xc[k];
    }

    // Read connectivity
    auto &ia= GetConnectivity();
    for (int k = 0; k<nelem*nvert_e; ++k)
    {
        ifs >> ia[k];
        ia[k] -= OFFSET;                // Matlab to C indexing
    }

    // additional read of boundary information (only start/end point)
    int nbedges;
    ifs >> nbedges;

/*
    bedges.resize(nbedges*2);
    for (int k = 0; k<nbedges*2; ++k)
    {
        ifs >> bedges[k];
        bedges[k] -= OFFSET;            // Matlab to C indexing
    }
*/   

    return;
}
