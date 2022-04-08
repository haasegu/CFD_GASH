// see:   http://llvm.org/docs/CodingStandards.html#include-style
#include "geom.h"
#include "utils.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
using namespace std;

Mesh::Mesh(int ndim, int nvert_e, int ndof_e)
    : _nelem(0), _nvert_e(nvert_e), _ndof_e(ndof_e), _nnode(0), _ndim(ndim), _ia(0), _xc(0)
{
}

Mesh::~Mesh()
{}

void Mesh::SetValues(std::vector<double> &v, const function<double(double,double)>& func) const
{
    int const nnode=Nnodes();              // number of vertices in mesh
    assert( nnode == static_cast<int>(v.size()) );
    for (int k=0; k<nnode; ++k)
    {
        v[k] = func( _xc[2*k], _xc[2*k+1] );
    }
}


void Mesh::Debug() const
{
    cout << "\n ############### Debug  M E S H  ###################\n";
    cout << "\n ...............    Coordinates       ...................\n";
    for (int k = 0; k < _nnode; ++k)
    {
        cout << k << " : ";
        for ( int i=0; i<_ndim; ++i )
        {
			cout << _xc[2 * k +i] << "  ";
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
    fout << Nnodes() << DELIMETER << Ndims() << DELIMETER << Nelems() << DELIMETER << NverticesElements() << endl;

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
    assert( Nelems()*NverticesElements() ==  static_cast<int>(_ia.size()));
    for (int k=0, kj=0; k<Nelems(); ++k)
    {
        for (int j=0; j<NverticesElements(); ++j, ++kj)
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





// #####################################################################
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



















