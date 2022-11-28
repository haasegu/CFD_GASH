// see:   http://llvm.org/docs/CodingStandards.html#include-style
#include "geom.h"
#include "utils.h"
#include "elements.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
using namespace std;

Mesh::Mesh(int ndim, int nvert_e, int ndof_e)
    : _nelem(0), _nvert_e(nvert_e), _ndof_e(ndof_e), _nnode(0), _ndim(ndim), _ia(0), _xc(0), _bedges(0), _dummy(0)
{
}

Mesh::~Mesh()
{}

void Mesh::SetValues(std::vector<double> &r, const function<double(double,double,double)>& func) const
{
    //int const nnode=Nnodes();              // number of vertices in mesh
    int const nnode=r.size();
    cout << "##  Mesh::SetValues : nnode " << nnode << "   size " << r.size() << endl;
    assert( nnode == static_cast<int>(r.size()) );
    for (int k=0; k<nnode; ++k)
    {
        r[k] = func( _xc[3*k], _xc[3*k+1] , _xc[3*k+2]);
    }
}

/*
void Mesh::SetValues1(std::vector<double> &u, const std::function<double(double, double, double)> &func1) const
{
    int const nnode=Nnodes();              // number of vertices in mesh
    cout << "##  Mesh::SetValues1 : nnode " << nnode << "   size " << u.size() << endl;
    //assert( nnode == static_cast<int>(u.size()) );
    for (int k=0; k<nnode; ++k)
    {
        u[k] = func1( _xc[3*k], _xc[3*k+1] , _xc[3*k+2]);
    }
}    

void Mesh::SetValues2(std::vector<double> &v, const std::function<double(double, double, double)> &func2) const
{
    int const nnode=Nnodes();              // number of vertices in mesh
    cout << "##  Mesh::SetValues2 : nnode " << nnode << "   size " << v.size() << endl;
    //assert( nnode == static_cast<int>(v.size()) );
    for (int k=0; k<nnode; ++k)
    {
        v[k] = func2( _xc[3*k], _xc[3*k+1] , _xc[3*k+2]);
    }
}        

void Mesh::SetValues3(std::vector<double> &w, const std::function<double(double, double, double)> &func3) const
{
    int const nnode=Nnodes();              // number of vertices in mesh
    cout << "##  Mesh::SetValues3 : nnode " << nnode << "   size " << w.size() << endl;
    //assert( nnode == static_cast<int>(w.size()) );
    for (int k=0; k<nnode; ++k)
    {
        w[k] = func3( _xc[3*k], _xc[3*k+1] , _xc[3*k+2]);
    }
}       
*/

void Mesh::Debug() const
{
    cout << "\n ############### Debug  M E S H  ###################\n";
    cout << "\n ...............    Coordinates       ...................\n";
    
    
    for (int k = 0; k < _nnode; ++k)
    {
        cout << k << " : ";
        for ( int i=0; i<_ndim; ++i )
        {
			cout << _xc[3 * k +i] << "  ";
		}
        cout << endl;
    }
    cout << "\n ...............    Elements        ...................\n";
    for (int k = 0; k < _nelem; ++k)
    {
        cout << k << " : ";
        for (int i = 0; i < _ndof_e; ++i )
            cout << _ia[_ndof_e * k + i] << "  ";
        cout << endl;
    }
    return;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||Periodic boundary conditions||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<double> Mesh::getPeriodicCoordsBox_xy(std::vector<int> &idx, double xl, double xh, double yl, double yh,double zl, double zh) const
{
	//const double BLA = -123456;
    std::vector<double> pCoords;// (Ndims()*idx.size(),BLA);
    std::vector<int> pidx;
    for(size_t k=0; k<idx.size(); k+=Ndims())
        {
		bool bp{true};
		double xp = _xc[idx[k]*Ndims()+0];
		double yp = _xc[idx[k]*Ndims()+1];
		double zp = _xc[idx[k]*Ndims()+2];
		if(equal(xp,xl)){
			xp=xh;}
		else if(equal(xp,xh)){
			xp=xl;}
		else if(equal(yp,yl)){
			yp=yh;}
		else if(equal(yp,yh)){
			yp=yl;}
		else {
//			 cout<<"Not found the periodic coordinate";
//			 assert(false);
                bp=false;
			 }
			 if (bp)
			 {
				 pCoords.push_back(xp);
				 pCoords.push_back(yp);
				 pCoords.push_back(zp);
				 pidx.push_back(idx[k]);
		     }
		}
		//auto ip=find(cbegin(pCoords),cend(pCoords),BLA);
		//assert(cend(pCoords)==ip);
		idx=pidx;
		//cout<<"pCoords = "<<pCoords<<endl;
		return pCoords;
}


std::vector<double> Mesh::getPeriodicCoordsBox_xz(std::vector<int> &idx, double xl, double xh, double yl, double yh,double zl, double zh) const
{
	//const double BLA = -123456;
    std::vector<double> pCoords;// (Ndims()*idx.size(),BLA);
    std::vector<int> pidx;
    for(size_t k=0; k<idx.size(); k+=Ndims())
        {
		bool bp{true};
		double xp = _xc[idx[k]*Ndims()+0];
		double yp = _xc[idx[k]*Ndims()+1];
		double zp = _xc[idx[k]*Ndims()+2];
		if(equal(xp,xl)){
			xp=xh;}
		else if(equal(xp,xh)){
			xp=xl;}
		else if(equal(zp,zl)){
			zp=zh;}
		else if(equal(zp,zh)){
			zp=zl;}
		else {
//			 cout<<"Not found the periodic coordinate";
//			 assert(false);
                bp=false;
			 }
			 if (bp)
			 {
				 pCoords.push_back(xp);
				 pCoords.push_back(yp);
				 pCoords.push_back(zp);
				 pidx.push_back(idx[k]);
		     }
		}
		//auto ip=find(cbegin(pCoords),cend(pCoords),BLA);
		//assert(cend(pCoords)==ip);
		idx=pidx;
		return pCoords;
}
std::vector<double> Mesh::getPeriodicCoordsBox_yz(std::vector<int> &idx, double xl, double xh, double yl, double yh,double zl, double zh) const
{
	//const double BLA = -123456;
    std::vector<double> pCoords;// (Ndims()*idx.size(),BLA);
    std::vector<int> pidx;
    for(size_t k=0; k<idx.size(); k+=Ndims())
        {
		bool bp{true};
		double xp = _xc[idx[k]*Ndims()+0];
		double yp = _xc[idx[k]*Ndims()+1];
		double zp = _xc[idx[k]*Ndims()+2];
		if(equal(yp,yl)){
			yp=yh;}
		else if(equal(yp,yh)){
			yp=yl;}
		else if(equal(zp,zl)){
			zp=zh;}
		else if(equal(zp,zh)){
			zp=zl;}
		else {
//			 cout<<"Not found the periodic coordinate";
//			 assert(false);
                bp=false;
			 }
			 if (bp)
			 {
				 pCoords.push_back(xp);
				 pCoords.push_back(yp);
				 pCoords.push_back(zp);
				 pidx.push_back(idx[k]);
		     }
		}
		//auto ip=find(cbegin(pCoords),cend(pCoords),BLA);
		//assert(cend(pCoords)==ip);
		idx=pidx;
		return pCoords;
}
double Mesh::getPeriodicValue(int k, std::vector<double> const &pCoords, std::vector<double> const &u) const
{
	    double xp = pCoords[Ndims()*k+0];
        double yp = pCoords[Ndims()*k+1];
        double zp = pCoords[Ndims()*k+2];
        
        int imin = -1; double dmin = 1e30;
        
        for(size_t i=0; i<_xc.size();i+=Ndims())
            {
			double disp = dist(xp,xp,zp,_xc[i+0],_xc[i+1],_xc[i+2]);
			if(dmin>disp)
			    {
					imin = i/Ndims();
					dmin = disp;
				}
			}
			return u[imin];
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||Periodic boundary conditions||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Mesh::Write_ascii_matlab(std::string const &fname, std::vector<double> const &v) const
{
    assert(Nnodes() ==  static_cast<int>(v.size()));  // fits vector length to mesh information?

    ofstream fout(fname);                             // open file ASCII mode
    if( !fout.is_open() )
    {
        cout << "\nFile " << fname << " has not been opened.\n\n" ;
        assert( fout.is_open() && "File not opened."  );
    }

    string const DELIMETER(" ");    // define the same delimeter as in matlab/ascii_read*.m
    int const    OFFSET(1);         // convert C-indexing to matlab

    // Write data: #nodes, #space dimensions, #elements, #vertices per element
    fout << Nnodes() << DELIMETER << Ndims() << DELIMETER << Nelems() << DELIMETER << NverticesElement() << endl;

    // Write cordinates: x_k, y_k   in seperate lines
    assert( Nnodes()*Ndims() ==  static_cast<int>(_xc.size()));
    for (int k=0, kj=0; k<Nnodes(); ++k)
    {
        for (int j=0; j<Ndims(); ++j, ++kj)
        {
            fout << _xc[kj] << DELIMETER;
        }
        fout << endl;
   }

    // Write connectivity: ia_k,0, ia_k,1 etc  in seperate lines
    assert( Nelems()*NverticesElement() ==  static_cast<int>(_ia.size()));
    for (int k=0, kj=0; k<Nelems(); ++k)
    {
        for (int j=0; j<NverticesElement(); ++j, ++kj)
        {
            fout << _ia[kj]+OFFSET << DELIMETER;       // C to matlab
        }
        fout << endl;
    }

    // Write vector
    for (int k=0; k<Nnodes(); ++k)
    {
        fout << v[k] << endl;
    }

    fout.close();
    return;
}


vector<vector<int>> Mesh::Node2NodeGraph_2() const
{
    return ::Node2NodeGraph(Nelems(),NdofsElement(),GetConnectivity());
}



vector<vector<int>> Node2NodeGraph(int const nelem, int const ndof_e,
                    vector<int> const &ia)
{
    assert(nelem*ndof_e==static_cast<int>(ia.size()));
    int const nnode = *max_element(cbegin(ia),cend(ia)) +1;
    vector<vector<int>> v2v(nnode, vector<int>(0));    // stores the vertex to vertex connections

    ////--------------
    vector<int> cnt(nnode,0);
    for (size_t i = 0; i < ia.size(); ++i)  ++cnt[ia[i]]; // determine number of entries per vertex
    for (size_t k = 0; k < v2v.size(); ++k)
    {
        v2v[k].resize(ndof_e * cnt[k]);              //    and allocate the memory for that vertex
        cnt[k] = 0;
    }
    ////--------------

    for (int e = 0; e < nelem; ++e)
    {
        int const basis = e * ndof_e;                  // start of vertex connectivity of element e
        for (int k = 0; k < ndof_e; ++k)
        {
            int const v = ia[basis + k];
            for (int l = 0; l < ndof_e; ++l)
            {
                v2v[v][cnt[v]] = ia[basis + l];
                ++cnt[v];
            }
        }
    }
    // finally  cnt[v]==v2v[v].size()  has to hold for all v!

    // guarantee unique, ascending sorted entries per vertex
    for (size_t v = 0; v < v2v.size(); ++v)
    {
        sort(v2v[v].begin(), v2v[v].end());
        auto ip = unique(v2v[v].begin(), v2v[v].end());
        v2v[v].erase(ip, v2v[v].end());
        //v2v[v].shrink_to_fit();       // automatically done when copied at return
    }

    return v2v;
}


void Mesh::Visualize(vector<double> const &v) const
{
    // define external command
    const string exec_m("matlab -nosplash < visualize_results.m");                 // Matlab
    //const string exec_m("octave --no-window-system --no-gui visualize_results.m"); // Octave
    //const string exec_m("flatpak run org.octave.Octave visualize_results.m");      // Octave (flatpak): desktop GH

    const string fname("uv.txt");
    Write_ascii_matlab(fname, v);

    int ierror=system(exec_m.c_str());                                   // call external command
    
    if (ierror!=0)
    {
        cout << endl << "Check path to Matlab/octave on your system" << endl;
    }
    cout << endl;
    return;
}

// subject to permutation:
//     re-sort:  _xc        
//               _xc[2*k_new], _xc[2*k_new+1]  with k_new = po2n[k] via old(_xc);
//     renumber: _ia, [_bedges, _edges ]
//                                 order ascending in each edge
//               old = _ia;
//               _ia[j] = p02n[old[j]]      j=0...3*ne-1
//
//     old2new = sort_indices(new_vertex_numbering)
void Mesh::PermuteVertices(std::vector<int> const& old2new)
{
    assert(Nnodes()==static_cast<int>(old2new.size()));

    permute_2(old2new, _xc);

    reNumberEntries(old2new, _ia);
  
    //reNumberEntries(old2new, _bedges);

    //reNumberEntries(old2new, _edges);
    //sortAscending_2(_edges);       // ascending order of vertices in edge
}
void Mesh::Write_ascii_paraview(std::string const &fname, std::vector<double> const &v) const
{
    assert(Nnodes() ==  static_cast<int>(v.size()));  // fits vector length to mesh information?

    ofstream fout(fname);                             // open file ASCII mode
    if ( !fout.is_open() )
    {
        cout << "\nFile " << fname << " has not been opened.\n\n" ;
        assert( fout.is_open() && "File not opened."  );
    }

    string const DELIMETER(" ");    // define the same delimiter as in matlab/ascii_read*.m
    //int const    OFFSET(o);         // C-indexing in output

    fout << "# vtk DataFile Version 2.0" << endl;
    fout << "HEAT EQUATION" << endl;
    fout << "ASCII" << endl;
    fout << "DATASET POLYDATA" << endl;
    fout << "POINTS "<< v.size()<<" float"<<endl;

    assert( Nnodes()*Ndims() ==  static_cast<int>(_xc.size()));
    for (int k = 0, kj = 0; k < Nnodes(); ++k)
    {
        for (int j = 0; j < Ndims(); ++j, ++kj)
        {
            fout << _xc[kj] << DELIMETER;
        }

        fout << v[k] << endl;
    }


    fout << "POLYGONS "<< Nelems() << ' ' << Nelems()*4 << endl;

    assert( Nelems()*NverticesElement() ==  static_cast<int>(_ia.size()));
    for (int k = 0, kj = 0; k < Nelems(); ++k)
    {
        fout << 3 << DELIMETER;          // triangular patches
        for (int j = 0; j < NverticesElement()-1; ++j, ++kj)
        {
            fout << _ia[kj] << DELIMETER;
        }

        fout << _ia[kj];
        kj=kj+1;

        if(k<Nelems()-1)
            {
                fout << endl;
            }
    }

    fout.close();
    return;
}
/*
 void Mesh::Write_ascii_paraview(std::string const &fname, std::vector<double> const &v) const
{
    //assert(3*Nnodes()+14 ==  static_cast<int>(v.size()));  // fits vector length to mesh information?

    ofstream fout(fname);                             // open file ASCII mode
    if ( !fout.is_open() )
    {
        cout << "\nFile " << fname << " has not been opened.\n\n" ;
        assert( fout.is_open() && "File not opened."  );
    }

    string const DELIMETER(" ");    // define the same delimiter as in matlab/ascii_read*.m
    //int const    OFFSET(o);         // C-indexing in output

    fout << "# vtk DataFile Version 3.0" << endl;
    fout << "HEAT EQUATION" << endl;
    fout << "ASCII" << endl;
    fout << "DATASET UNSTRUCTURED_GRID" << endl;
    fout << "POINTS "<< v.size()<<" float"<<endl;

    assert( Nnodes()*Ndims() ==  static_cast<int>(_xc.size()));
    for (int k = 0, kj = 0; k < Nnodes(); ++k)
    {
        for (int j = 0; j < Ndims(); ++j, ++kj)
        {
            fout << _xc[kj] << DELIMETER;
        }

        fout << endl;
    }


    fout << "CELLS "<< Nelems() << ' ' << Nelems()*4 << endl;

    assert( Nelems()*NverticesElement() ==  static_cast<int>(_ia.size()));
    for (int k = 0, kj = 0; k < Nelems(); ++k)
    {
        fout << 4 << DELIMETER;          // triangular patches
        for (int j = 0; j < NverticesElement()-1; ++j, ++kj)
        {
            fout << _ia[kj] << DELIMETER;
        }

        fout << _ia[kj];
        kj=kj+1;

        if(k<Nelems()-1)
            {
                fout << endl;
            }
    }
    fout << endl;
    fout << "CELL_TYPES "<< Nelems() <<endl;
    for (int k=0; k<Nelems();++k){
		
		fout << 10 << endl;
		
		}
	
    fout << "CELL_DATA "<< Nelems() <<endl;
    fout << "POINT_DATA "<< v.size() <<endl;
    fout << "SCALARS temp float" <<endl;
    fout << "LOOKUP_TABLE default" <<endl;

    assert( Nnodes()*Ndims() ==  static_cast<int>(_xc.size()));
    for (int k = 0; k < Nnodes(); ++k)
    {
        fout << v[k]<< endl;
    }
  

    fout.close();
    return;
}
*/
void Mesh::Visualize_paraview(vector<double> const &v) const
{
    //const string exec_m("open -a paraview");                 // paraview
    const string exec_m("paraview");                 // paraview
   
    const string fname("uv.vtk");
    Write_ascii_paraview(fname, v);

    int ierror = system(exec_m.c_str());                                 // call external command

    if (ierror != 0)
    {
        cout << endl << "Check path to paraview on your system" << endl;
    }
    cout << endl;
    return;
} 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||Boundary nodes/conditions||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<int> Mesh::Index_DirichletNodes() const
{
    //vector<int> idx(_bedges);                              // copy

    //sort(idx.begin(), idx.end());                         // sort
    //idx.erase( unique(idx.begin(), idx.end()), idx.end() ); // remove duplicate data

    //return idx;
    vector<int> idx;
    cout << "\n WARNING:  Mesh::Index_DirichletNodes()   N O T   implemented\n";
    return idx;
} 

vector<int> Mesh::Index_DirichletNodes_Box
    (double xl, double xh, double yl, double yh) const
{
    assert(2==Ndims());        // not in 3D currently
	auto x=GetCoords();
	vector<int> idx;
	for (int k=0; k<Nnodes()*Ndims(); k+=2)
	{
		const double xk(x[k]), yk(x[k+1]);
		if (equal(xk,xl) || equal(xk,xh) || equal(yk,yl) || equal(yk,yh))
		{
			idx.push_back(k/2);
		}
	}
    
    sort(idx.begin(), idx.end());                           // sort
    idx.erase( unique(idx.begin(), idx.end()), idx.end() ); // remove duplicate data    
	return idx;
}

vector<int> Mesh::Index_DirichletNodes_Box
    (double xl, double xh, double yl, double yh,double zl, double zh) const
{
    assert(3==Ndims());        // not in 3D currently
	auto x=GetCoords();
	vector<int> idx;
	for (int k=0; k<Nnodes()*Ndims(); k+=3)
	{
		const double xk(x[k]), yk(x[k+1]), zk(x[k+2]);
		if (equal(xk,xl) || equal(xk,xh) || equal(yk,yl) || equal(yk,yh) ||  equal(zk,zl) || equal(zk,zh) )
		{
			idx.push_back(k/3);
		}
	}
    
    sort(idx.begin(), idx.end());                           // sort
    idx.erase( unique(idx.begin(), idx.end()), idx.end() ); // remove duplicate data    
	return idx;
}


Mesh::Mesh(std::string const &fname)
    //: Mesh(2, 3, 3, 3) // two dimensions, 3 vertices, 3 dofs, 3 edges per element
    : Mesh()
{
    ReadVertexBasedMesh(fname);
    //Debug(); int ijk; cin >>ijk;
    
    //liftToQuadratic();
    //Debug();
    
    //DeriveEdgeFromVertexBased();        // Generate also the edge based information
}

void Mesh::ReadVertexBasedMesh(std::string const &fname)
{
    ifstream ifs(fname);
    if (!(ifs.is_open() && ifs.good()))
    {
        cerr << "Mesh::ReadVertexBasedMesh: Error cannot open file " << fname << endl;
        assert(ifs.is_open());
    }

    int const OFFSET(1);             // Matlab to C indexing
    cout << "ASCI file  " << fname << "  opened" << endl;

    // Read some mesh constants
    int nnode, ndim, nelem, nvert_e;
    ifs >> nnode >> ndim >> nelem >> nvert_e;
    cout << nnode << "  " << ndim << "  " << nelem << "  " << nvert_e << endl;
    // accept only    triangles (2D)  or  tetrahedrons (3D)
    assert((ndim == 2 && nvert_e == 3)||(ndim == 3 && nvert_e == 4));
    
    // set member
    _ndim    = ndim;
    _nvert_e = nvert_e;
    _ndof_e  = _nvert_e;
    //_nedge_e = (2==_ndim)? 3:6;       // Generate also the edge based information

    // Allocate memory
    Resize_Coords(nnode, ndim);                 // coordinates in 2D [nnode][ndim]
    Resize_Connectivity(nelem, nvert_e);        // connectivity matrix [nelem][nvert_e]

    // Read coordinates
    auto &xc = GetCoords();
    for (int k = 0; k < nnode * ndim; ++k)
    {
        ifs >> xc[k];
    }

    // Read connectivity
    auto &ia = GetConnectivity();
    for (int k = 0; k < nelem * nvert_e; ++k)
    {
        ifs >> ia[k];
        ia[k] -= OFFSET;                // Matlab to C indexing
    }

    if (2==ndim)
    {
        // additional read of (Dirichlet) boundary information (only start/end point)
        int nbedges;
        ifs >> nbedges;

        _bedges.resize(nbedges * 2);
        for (int k = 0; k < nbedges * 2; ++k)
        {
            ifs >> _bedges[k];
            _bedges[k] -= OFFSET;            // Matlab to C indexing
        }
    }
    else
    {
        // ToDo: add boundary information to 3D mesh
        cout << std::endl << "NO boundary information available for 3D mesh" << endl;
    }
    return;
}


void Mesh::liftToQuadratic()
{
    cout << "##  Mesh::liftToQuadratic  ##" << endl;
    int const nelem   = Nelems();           // number of elements remains unchanged
    int const nnodes1 = Nnodes();           // number of P1-vertices 
    int const nvert_1 = NverticesElement(); // #vertices per P1-element
    assert(NverticesElement()==NdofsElement());
    
    vector<double> const xc_p1 = GetCoords();     // save P1 coordinates
    vector<int> const ia_p1 = GetConnectivity();  // save P1 connevctivity
    // check dimensions in P1
    assert( nnodes1*Ndims()==static_cast<int>(xc_p1.size()) );
    assert( nelem*nvert_1==static_cast<int>(ia_p1.size()) );
    
    bool lifting_possible{false};
    // P1 --> P2 DOFSs per element
    if (2==Ndims())
    {
        lifting_possible = (3==nvert_1);
        if (lifting_possible)
        { 
            SetNverticesElement(6);
            SetNdofsElement(NverticesElement());
        }
//        assert(false);
    }
    else if(3==Ndims())
    {
        lifting_possible = (4==nvert_1);
        if (lifting_possible)
        { 
            SetNverticesElement(10);
            SetNdofsElement(NverticesElement());
        }
    }
    else
    {
        cout << "Mesh::liftToQuadratic(): Wrong space dimension :" << Ndims() << endl;
        assert(false);
    }
    
    if (!lifting_possible) 
    {
        cout << "Mesh::liftToQuadratic(): Mesh elements must be linear." << endl;
        return;                                   // Exit function and do nothing.
    }
    
    int const nvert_2 = NverticesElement();       // #vertices per P2-element
    cout << "nvert_2: " << nvert_2 << endl;

    // P2 connectivity: memory allocation and partial initialization with P1 data
    vector<int>  & ia_p2 = GetConnectivity();     // P2 connectivity
    ia_p2.resize(nelem*nvert_2,-1);
    // Copy P1 connectivity into P2 
    for (int ke=0; ke<nelem; ++ke)
    {
        int const idx1=ke*nvert_1;
        int const idx2=ke*nvert_2;

        for (int d=0; d<nvert_1; ++d)
        {
            ia_p2.at(idx2+d) = ia_p1.at(idx1+d);
        }
    }
    //cout << "SALMAN: " << endl;
    // P2 coordinates: max. memory reservation and partial initialization with P1 data
    vector<double>  & xc_p2 = GetCoords();
    // reserve max. memory, append new vertices with push_back(),  call shrink_to_fit() finally.
    xc_p2.reserve(nnodes1+(nvert_2-nvert_1)*nelem);
    xc_p2.resize(nnodes1*Ndims(),-12345);
    copy(cbegin(xc_p1),cend(xc_p1),begin(xc_p2));
    //cout << "SALMAN: " << endl;

    int off_2d = 3-Ndims();
    for (int ke=0; ke<nelem; ++ke)
    {
        int const idx2=ke*nvert_2;                // Element starts
        int const v0 = ia_p2.at(idx2+0);          // vertices of P1
        int const v1 = ia_p2.at(idx2+1);
        int const v2 = ia_p2.at(idx2+2);
        //cout << v0 << "  " << v1 << "  " << v2 << "  ::  " << idx2 << "  " << ia_p2.size() <<"\n";
        
        // GH: not correct in 2D!
        ia_p2.at(idx2+4-off_2d) = appendMidpoint(v0,v1,xc_p2);
        ia_p2.at(idx2+5-off_2d) = appendMidpoint(v1,v2,xc_p2);
        ia_p2.at(idx2+6-off_2d) = appendMidpoint(v2,v0,xc_p2);
        //cout <<" v0: "<<ia_p2.at(idx2+0)<<" v1: "<<ia_p2.at(idx2+1)<<" v2: "<<ia_p2.at(idx2+2)<<" v3:  "<<ia_p2.at(idx2+4-off_2d)<<" v4:  "<<ia_p2.at(idx2+5-off_2d)<<" v5:  "<<ia_p2.at(idx2+6-off_2d)<<endl;
        if (3==Ndims())
        {
            int const v3 = ia_p2.at(idx2+3);      // forth vertex of P1 in 3D
            ia_p2.at(idx2+7) = appendMidpoint(v0,v3,xc_p2);
            ia_p2.at(idx2+8) = appendMidpoint(v1,v3,xc_p2);
            ia_p2.at(idx2+9) = appendMidpoint(v2,v3,xc_p2);
        }
    }
    
    //cout << "SALMAN: " << endl;
    xc_p2.shrink_to_fit();
    SetNnode(xc_p2.size()/Ndims()); 
    
    cout << _nnode << "  " << _ndim << "  " << _nelem << "  " << _nvert_e << "  " << _ndof_e << endl;

}


// #######################################################################################################################################################################################################
Mesh_2d_3_square::Mesh_2d_3_square(int nx, int ny, int myid, int procx, int procy)
    : Mesh(2,3,3),  // two dimensions, 3 vertices, 3 dofs
      _myid(myid), _procx(procx), _procy(procy), _neigh{{-1,-1,-1,-1}}, _color(0),
      _xl(0.0), _xr(1.0), _yb(0.0), _yt(1.0), _nx(nx), _ny(ny)
{
    //void IniGeom(int const myid, int const procx, int const procy, int neigh[], int &color)
    int const ky = _myid / _procx;
    int const kx = _myid % _procy;	//    MOD(myid,procx)
    // Determine the neighbors of domain/rank myid
    _neigh[0] = (ky == 0)       ?  -1 : _myid - _procx;    //   South
    _neigh[1] = (kx == _procx - 1) ?  -1 : _myid + 1;      //   East
    _neigh[2] = (ky == _procy - 1) ?  -1 : _myid + _procx;  //   North
    _neigh[3] = (kx == 0)       ?  -1 : _myid - 1;        //   West

    _color = (kx + ky) & 1 ;

    // subdomain is part of unit square
    double const hx = 1. / _procx;
    double const hy = 1. / _procy;
    _xl = kx * hx;                      //  left
    _xr = (kx + 1) * hx;                //  right
    _yb = ky * hy;                      //  bottom
    _yt = (ky + 1) * hy;                //  top

    // Calculate coordinates
    int const nnode = (_nx + 1) * (_ny + 1); // number of nodes
    Resize_Coords(nnode, 2);                 // coordinates in 2D [nnode][ndim]
    GetCoordsInRectangle(_nx, _ny, _xl, _xr, _yb, _yt, GetCoords().data());

    // Calculate element connectivity (linear triangles)
    int const nelem = 2 * _nx * _ny;         // number of elements
    Resize_Connectivity(nelem, 3);           // connectivity matrix [nelem][3]
    GetConnectivityInRectangle(_nx, _ny, GetConnectivity().data());
    
        
    //// GH
    //vector<int> perm(Nnodes());
    //iota(rbegin(perm),rend(perm),0);
    ////random_shuffle(begin(perm),end(perm));
    //PermuteVertices(perm);

    return;
}


void Mesh_2d_3_square::SetU(std::vector<double> &u) const
{
    int dx    = _nx + 1;
    for (int j = 0; j <= _ny; ++j)
    {
        int k = j * dx;
        for (int i = 0; i <= _nx; ++i, ++k)
        {
            u[k] = 0.0;
        }
    }

}

void Mesh_2d_3_square::SetF(std::vector<double> &f) const
{
    int dx    = _nx + 1;
    for (int j = 0; j <= _ny; ++j)
    {
        int k = j * dx;
        for (int i = 0; i <= _nx; ++i, ++k)
        {
            f[k] = 1.0;
        }
    }

}


std::vector<int> Mesh_2d_3_square::Index_DirichletNodes() const
{
    int const dx = 1,
              dy = _nx + 1,
              bl  = 0,
              br  = _nx,
              tl  = _ny * (_nx + 1),
              tr  = (_ny + 1) * (_nx + 1) - 1;
    int const start[4] = { bl, br, tl, bl},
                end[4] = { br, tr, tr, tl},
               step[4] = { dx, dy, dx, dy};

    vector<int> idx(0);
    for (int j = 0; j < 4; j++)
    {
        if (_neigh[j] < 0)
        {
            for (int i = start[j]; i <= end[j]; i += step[j])
            {
                idx.push_back(i);        // node i is Dirichlet node
            }
        }
    }
    //    remove multiple elements
    sort(idx.begin(),idx.end());                           // sort
    idx.erase( unique(idx.begin(),idx.end()), idx.end() ); // remove duplicate data

    return idx;
}


void Mesh_2d_3_square::SaveVectorP(std::string const &name, vector<double> const &u) const
{
//  construct the file name for subdomain myid
    const string tmp( std::to_string(_myid / 100) + to_string((_myid % 100) / 10) + to_string(_myid % 10) );

    const string namep = name + "." + tmp;
    ofstream ff(namep.c_str());
    ff.precision(6);
    ff.setf(ios::fixed, ios::floatfield);

    // assumes tensor product grid in unit square; rowise numbered (as generated in class constructor)
    // output is provided for tensor product grid visualization ( similar to Matlab-surf() )
    auto const &xc = GetCoords();
    int k = 0;
    for (int j = 0; j <= _ny; ++j)
    {
        for (int i = 0; i <= _nx; ++i, ++k)
            ff << xc[2 * k + 0] << "   " << xc[2 * k + 1] << "   " << u[k] << endl;
        ff << endl;
    }

    ff.close();
    return;
}

void Mesh_2d_3_square::GetCoordsInRectangle(int const nx, int const ny,
                          double const xl, double const xr, double const yb, double const yt,
                          double xc[])
{
    const double hx = (xr - xl) / nx,
                 hy = (yt - yb) / ny;

    int k = 0;
    for (int j = 0; j <= ny; ++j)
    {
        const double y0 = yb + j * hy;
        for (int i = 0; i <= nx; ++i, k += 2)
        {
            xc[k  ] = xl + i * hx;
            xc[k + 1] = y0;
        }
    }

    return;
}

void Mesh_2d_3_square::GetConnectivityInRectangle(int const nx, int const ny, int ia[])
{
    const int dx = nx + 1;
    int k  = 0;
    int l  = 0;
    for (int j = 0; j < ny; ++j, ++k)
    {
        for (int i = 0; i < nx; ++i, ++k)
        {
            ia[l  ] = k;
            ia[l + 1] = k + 1;
            ia[l + 2] = k + dx + 1;
            l += 3;
            ia[l  ] = k;
            ia[l + 1] = k + dx;
            ia[l + 2] = k + dx + 1;
            l += 3;
        }
    }
    return;
}

// #################### still some old code (--> MPI) ############################


//  Copies the values of w corresponding to the boundary
//  South (ib==1), East (ib==2), North (ib==3), West (ib==4)

void GetBound(int const ib, int const nx, int const ny, double const w[], double s[])
{
    const int //dx = 1,
    dy = nx + 1,
    bl  = 0,
    br  = nx,
    tl  = ny * (nx + 1),
    tr  = (ny + 1) * (nx + 1) - 1;
    switch (ib)
    {
        case 1:
        {
            for (int i = bl, j = 0; i <= br; ++i, ++j)
                s[j] = w[i];
            break;
        }
        case 3:
        {
            for (int i = tl, j = 0; i <= tr; ++i, ++j)
                s[j] = w[i];
            break;
        }
        case 4:
        {
            for (int i = bl, j = 0; i <= tl; i += dy, ++j)
                s[j] = w[i];
            break;
        }
        case 2:
        {
            for (int i = br, j = 0; i <= tr; i += dy, ++j)
                s[j] = w[i];
            break;
        }
        default:
        {
            cout << endl << "Wrong parameter  ib in " << __FILE__ << ":" << __LINE__ << endl;
        }
    }
    return;
}

// ----------------------------------------------------------------------------------------------------------
// Computes w:  = w + s at nodes on  the boundary
// South (ib == 1), East (ib == 2), North (ib == 3), West (ib == 4)

void AddBound(int const ib, int const nx, int const ny, double w[], double const s[])
{
    int const dy = nx + 1,
              bl  = 0,
              br  = nx,
              tl  = ny * (nx + 1),
              tr  = (ny + 1) * (nx + 1) - 1;
    switch (ib)
    {
        case 1:
        {
            for (int i = bl, j = 0; i <= br; ++i, ++j)
                w[i] += s[j];
            break;
        }
        case 3:
        {
            for (int i = tl, j = 0; i <= tr; ++i, ++j)
                w[i] += s[j];
            break;
        }
        case 4:
        {
            for (int i = bl, j = 0; i <= tl; i += dy, ++j)
                w[i] += s[j];
            break;
        }
        case 2:
        {
            for (int i = br, j = 0; i <= tr; i += dy, ++j)
                w[i] += s[j];
            break;
        }
        default:
        {
            cout << endl << "Wrong parameter  ib in " << __FILE__ << ":" << __LINE__ << endl;
        }
    }
    return;
}


// ####################################################################

Mesh_2d_3_matlab::Mesh_2d_3_matlab(string const &fname)
 : Mesh(2,3,3),    // two dimensions, 3 vertices, 3 dofs
   bedges(0)
{
    ifstream ifs(fname);
    if(!(ifs.is_open() && ifs.good()))
    {
        cerr << "Mesh_2d_3_matlab: Error cannot open file " << fname << endl;
        assert(ifs.is_open());
    }

    int const OFFSET(1);             // Matlab to C indexing
    cout << "ASCI file  " << fname << "  opened" << endl;

    // Read some mesh constants
    int nnode, ndim, nelem, nvert_e;
    ifs >> nnode >> ndim >> nelem >> nvert_e;
    cout << nnode << "  " << ndim << "  " << nelem << "  " << nvert_e <<endl;
    assert(ndim==2 && nvert_e==3);

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

    bedges.resize(nbedges*2);
    for (int k = 0; k<nbedges*2; ++k)
    {
        ifs >> bedges[k];
        bedges[k] -= OFFSET;            // Matlab to C indexing
    }
    
    // GH
    cout << " REORDER  VERTICES " << endl;
    vector<int> perm(Nnodes());
    iota(rbegin(perm),rend(perm),0);
    random_shuffle(begin(perm),end(perm));
    PermuteVertices(perm);


    return;
}

// binary
//{
    //ifstream ifs(fname, ios_base::in | ios_base::binary);
    //if(!(ifs.is_open() && ifs.good()))
    //{
        //cerr << "ReadBinMatrix: Error cannot open file " << file << endl;
        //assert(ifs.is_open());
    //}
    //cout << "ReadBinMatrix: file opened" << file << endl;


//}

// binaryIO.cpp
//void read_binMatrix(const string& file, vector<int> &cnt, vector<int> &col, vector<double> &ele)
//{

    //ifstream ifs(file, ios_base::in | ios_base::binary);

    //if(!(ifs.is_open() && ifs.good()))
    //{
        //cerr << "ReadBinMatrix: Error cannot open file " << file << endl;
        //assert(ifs.is_open());
    //}
    //cout << "ReadBinMatrix: Opened file " << file << endl;

    //int _size;

    //ifs.read(reinterpret_cast<char*>(&_size), sizeof(int));   // old: ifs.read((char*)&_size, sizeof(int));
    //cnt.resize(_size);
    //cout << "ReadBinMatrix: cnt size: " << _size << endl;

    //ifs.read((char*)&_size, sizeof(int));
    //col.resize(_size);
    //cout << "ReadBinMatrix: col size: " << _size << endl;

    //ifs.read((char*)&_size, sizeof(int));
    //ele.resize(_size);
    //cout << "ReadBinMatrix: ele size: " << _size << endl;


    //ifs.read((char*)cnt.data(), cnt.size() * sizeof(int));
    //ifs.read((char*)col.data(), col.size() * sizeof(int));
    //ifs.read((char*)ele.data(), ele.size() * sizeof(double));

    //ifs.close();
    //cout << "ReadBinMatrix: Finished reading matrix.." << endl;

//}

void Mesh_2d_3_matlab::PermuteVertices(std::vector<int> const& old2new)
{
    Mesh::PermuteVertices(old2new);
    reNumberEntries(old2new, bedges);

    //reNumberEntries(old2new, _edges);
    //sortAscending_2(_edges);       // ascending order of vertices in edge
}

std::vector<int> Mesh_2d_3_matlab::Index_DirichletNodes() const
{
    vector<int> idx(bedges);                              // copy

    sort(idx.begin(),idx.end());                          // sort
    idx.erase( unique(idx.begin(),idx.end()), idx.end() );// remove duplicate data

    return idx;
}



















