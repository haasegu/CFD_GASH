
#include "geom3.h"
#include "vdop.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
using namespace std;





// ####################################################################

Mesh_3d_4_matlab::Mesh_3d_4_matlab(string const &fname)
 : Mesh(3,4,4)    // three dimensions, 4 vertices, 4 dofs
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

std::vector<int>  Mesh_3d_4_matlab::Index_DirichletNodes_Box1(double xl, double xh, double yl, double yh,double zl, double zh )
{
	auto x=GetCoords();
	std::vector<int> iDir;
	for (int k=0; k<Nnodes()*Ndims(); k+=3)
	{
		const double xk(x[k]), yk(x[k+1]), zk(x[k+2]);
		if (equal(xk,xl) || equal(xk,xh) || equal(yk,yl) || equal(yk,yh) ||  equal(zk,zl) || equal(zk,zh) )
		{
			iDir.push_back(k/3);
		}
	}
	return iDir;
	
}


// ####################################################################

Mesh_3d_P1_matlab::Mesh_3d_P1_matlab(string const &fname)
 : Mesh(3,4,4)        // three dimensions, 4 vertices, 4 dofs
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
std::vector<int>  Mesh_3d_P1_matlab::Index_DirichletNodes_Box(double xl, double xh, double yl, double yh,double zl, double zh )
{
	auto x=GetCoords();
	std::vector<int> iDir;
	for (int k=0; k<Nnodes()*Ndims(); k+=3)
	{
		const double xk(x[k]), yk(x[k+1]), zk(x[k+2]);
		if (equal(xk,xl) || equal(xk,xh) || equal(yk,yl) || equal(yk,yh) ||  equal(zk,zl) || equal(zk,zh) )
		{
			iDir.push_back(k/3);
		}
	}
	return iDir;
	
}
*/
// ####################################################################

Mesh_3d_P2_matlab::Mesh_3d_P2_matlab(string const &fname)
 : Mesh(3,10,10)       // three dimensions, 10 vertices, 10 dofs
 //: Mesh_3d_P1_matlab(fname)
   //bedges(0)
{
    Mesh_3d_P1_matlab const p1_mesh(fname);        // read geometry of P1 mesh
    deriveMeshFromP1(p1_mesh);               // --> P2 mesh

    return;
}


void Mesh_3d_P2_matlab::deriveMeshFromP1(Mesh_3d_P1_matlab const &p1)
{
    Mesh_3d_P2_matlab &p2 = *this;
    
    //p1.Debug();

    int nelem = p1.Nelems();
    p2.SetNelem(nelem);                // number of elements
    int nnodes1 = p1.Nnodes();
    //p2.SetNnode(nnodes1);                // number of vertices   
    
     // check dimensions in P1
    vector<double> const & xc_p1 = p1.GetCoords();
    assert( nnodes1*3==static_cast<int>(xc_p1.size()) );
    vector<int> const & ia_p1 = p1.GetConnectivity();
    assert( nelem*4==static_cast<int>(ia_p1.size()) );

    // P2 connectivity: memory allocation and partial initialization with P1 data
    vector<int>       & ia_p2 = p2.GetConnectivity();
    ia_p2.resize(nelem*p2.NdofsElement(),-1);
    // Copy P1 connectivity into P2 
    for (int ke=0; ke<nelem; ++ke)
    {
        int idx1=ke*p1.NdofsElement();
        int idx2=ke*p2.NdofsElement();

        for (int d=0; d<p1.NdofsElement(); ++d)
        {
            ia_p2.at(idx2+d) = ia_p1.at(idx1+d);
        }
    }
    
    // p2 coordinates: max. memory reservation and partial initialization with P1 data
    vector<double>  & xc_p2 = p2.GetCoords();
    // reserve max. memory, append new vertices with push_back(),  call shrink_to_fit() finally.
    xc_p2.reserve(nnodes1+(p2.NverticesElement()-p1.NverticesElement())*nelem);
    xc_p2.resize(nnodes1*3,-12345);
    copy(cbegin(xc_p1),cend(xc_p1),begin(xc_p2));
    
    //int nnodes2=nnodes1-1;                          // we are going to generate new vertices
    for (int ke=0; ke<nelem; ++ke)
    {
        int const idx2=ke*p2.NdofsElement();      // Element starts
        int const v0 = ia_p2.at(idx2+0);
        int const v1 = ia_p2.at(idx2+1);
        int const v2 = ia_p2.at(idx2+2);
        int const v3 = ia_p2.at(idx2+3);          // vertices of P1
        int idxv;
        
        idxv = appendMidpoint(v0,v1,xc_p2);
        ia_p2.at(idx2+4) = idxv;
        idxv = appendMidpoint(v1,v2,xc_p2);
        ia_p2.at(idx2+5) = idxv;
        idxv = appendMidpoint(v2,v0,xc_p2);
        ia_p2.at(idx2+6) = idxv;
        
        idxv = appendMidpoint(v0,v3,xc_p2);
        ia_p2.at(idx2+7) = idxv;
        idxv = appendMidpoint(v1,v3,xc_p2);
        ia_p2.at(idx2+8) = idxv;
        idxv = appendMidpoint(v2,v3,xc_p2);
        ia_p2.at(idx2+9) = idxv;
           
    }
    
    xc_p2.shrink_to_fit();
    
    p2.SetNnode(xc_p2.size()/3); p2.Debug();
}

int appendMidpoint(int v1, int v2, vector<double> &xc, int ndim)
{
//	cout << "appendMidpoint\n";
    assert(2==ndim);            // works also for 2D
    int const i1{v1*ndim};
    int const i2{v2*ndim};
    vector<double> const xm{(xc.at(i1+0)+xc.at(i2+0))/2, (xc.at(i1+1)+xc.at(i2+1))/2};
    int idx_vertex=getVertexIndex(xm, xc);
    if (0>idx_vertex)
    {
        for (int d=0; d<ndim; ++d)
        {
            xc.push_back(xm[d]);
        }
        idx_vertex = static_cast<int>(xc.size())/ndim-1;
    }
    
    return idx_vertex;
}

int getVertexIndex(vector<double> const &xm, vector<double> const &xc, int ndim)
{
    assert(ndim==static_cast<int>(xm.size()));            // works also for 2D
//    cout << "getVertexIndex" << endl;
    auto xcStart{cbegin(xc)};
    int idx=-1;
    size_t k=0;
    while (idx<0 && k<xc.size())
    {
        //if ( equal( cbegin(xm),cend(xm), xcStart+k ) )
        if ( equal( cbegin(xm),cend(xm), xcStart+k, [](double a, double b){return equal(a,b);} ) )
        {
            idx = static_cast<int>(k)/ndim;
        }
        k+=ndim;
    }
    return idx;
}



std::vector<int>  Mesh_3d_P2_matlab::Index_DirichletNodes_Box1(double xl, double xh, double yl, double yh,double zl, double zh )
{
	auto x=GetCoords();
	std::vector<int> iDir;
	for (int k=0; k<Nnodes()*Ndims(); k+=3)
	{
		const double xk(x[k]), yk(x[k+1]), zk(x[k+2]);
		if (equal(xk,xl) || equal(xk,xh) || equal(yk,yl) || equal(yk,yh) ||  equal(zk,zl) || equal(zk,zh) )
		{
			iDir.push_back(k/3);
		}
	}
	return iDir;
	
}



