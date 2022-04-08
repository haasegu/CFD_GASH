
#include "geom3.h"
#include "vdop.h"

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




/*
gMesh_Hierarchy::gMesh_Hierarchy(Mesh const &cmesh, int const nlevel)
    : _gmesh(max(1, nlevel))
{
    _gmesh[0] = make_shared<Mesh>(cmesh);
    for (int lev = 1; lev < nlevel; ++lev)
    {
        _gmesh.at(lev) = make_shared<RefinedMesh>( *_gmesh.at(lev - 1) );
        //auto vv=_gmesh[lev]->GetFathersOfVertices();
        //cout << " :: "<< vv.size() <<endl;
    }
    for (size_t lev = 0; lev < _gmesh.size(); ++lev)
    {
        _gmesh[lev]->Del_EdgeConnectivity();
    }
}
*/

std::vector<int>  Mesh_3d_4_matlab::Index_DirichletNodes_Box(double xl, double xh, double yl, double yh,double zl, double zh )
{
	auto x=GetCoords();
	std::vector<int> iDir;
	for (int k=0; k<Nnodes()*Ndims(); k+=3)
	{
		const double xk(x[k]), yk(x[k+1]), zk(x[k+2]);
		if (xk==xl || xk==xh || yk==yl || yk==yh ||  zk==zl || zk==zh )
		{
			iDir.push_back(k/3);
		}
	}
	return iDir;
	
}
