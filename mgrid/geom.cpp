// see:   http://llvm.org/docs/CodingStandards.html#include-style
#include "vdop.h"
#include "cuthill_mckee_ordering.h"
#include "geom.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <ctime>                  // contains clock()
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

using namespace std;

Mesh::Mesh(int ndim, int nvert_e, int ndof_e)
    : _nelem(0), _nvert_e(nvert_e), _ndof_e(ndof_e), _nnode(0), _ndim(ndim), _ia(0), _xc(0),
      _bedges(0),
      _nedge(0), _edges(0), _ea(), _ebedges(),
      _dummy(0)
{
}

Mesh::~Mesh()
{}

void Mesh::SetValues(std::vector<double> &v, const function<double(double, double)> &func) const
{
    int const nnode = Nnodes();            // number of vertices in mesh
    assert( nnode == static_cast<int>(v.size()) );
    for (int k = 0; k < nnode; ++k)
    {
        v[k] = func( _xc[2 * k], _xc[2 * k + 1] );
    }
}


void Mesh::Debug() const
{
    cout << "\n ############### Debug  M E S H  ###################\n";
    cout << "\n ...............    Coordinates       ...................\n";
    for (int k = 0; k < _nnode; ++k)
    {
        cout << k << " : " << _xc[2 * k] << "  " << _xc[2 * k + 1] << endl;
    }
    cout << "\n ...............    Elements        ...................\n";
    for (int k = 0; k < _nelem; ++k)
    {
        cout << k << " : ";
        for (int i = 0; i < _ndof_e; ++i )
            cout << _ia[_ndof_e * k + i] << "  ";
        cout << endl;
    }
    cout << "\n ...............    Boundary (vertices)    .................\n";
    cout << " _bedges : " << _bedges << endl;
    return;
}

void Mesh::DebugEdgeBased() const
{
    cout << "\n ############### Debug  M E S H  (edge based) ###################\n";
    cout << "\n ...............    Coordinates       ...................\n";
    for (int k = 0; k < _nnode; ++k)
    {
        cout << k << " : " << _xc[2 * k] << "  " << _xc[2 * k + 1] << endl;
    }

    cout << "\n ...............    edges        ...................\n";
    for (int k = 0; k < _nedge; ++k)
    {
        cout << k << " : ";
        for (int i = 0; i < 2; ++i )
            cout << _edges[2 * k + i] << "  ";
        cout << endl;
    }

    cout << "\n ...............    Elements (edges)    .................\n";
    assert(_nedge_e * _nelem == static_cast<int>(_ea.size()) );
    for (int k = 0; k < _nelem; ++k)
    {
        cout << k << " : ";
        for (int i = 0; i < _nedge_e; ++i )
            cout << _ea[_nedge_e * k + i] << "  ";
        cout << endl;
    }
    cout << "\n ...............    Boundary (edges)    .................\n";
    cout << " _ebedges : " << _ebedges << endl;

    return;
}

void Mesh::Write_ascii_matlab(std::string const &fname, std::vector<double> const &v) const
{
    assert(Nnodes() ==  static_cast<int>(v.size()));  // fits vector length to mesh information?

    ofstream fout(fname);                             // open file ASCII mode
    if ( !fout.is_open() )
    {
        cout << "\nFile " << fname << " has not been opened.\n\n" ;
        assert( fout.is_open() && "File not opened."  );
    }

    string const DELIMETER(" ");    // define the same delimiter as in matlab/ascii_read*.m
    int const    OFFSET(1);         // convert C-indexing to matlab

    // Write data: #nodes, #space dimensions, #elements, #vertices per element
    fout << Nnodes() << DELIMETER << Ndims() << DELIMETER << Nelems() << DELIMETER << NverticesElements() << endl;

    // Write coordinates: x_k, y_k   in separate lines
    assert( Nnodes()*Ndims() ==  static_cast<int>(_xc.size()));
    for (int k = 0, kj = 0; k < Nnodes(); ++k)
    {
        for (int j = 0; j < Ndims(); ++j, ++kj)
        {
            fout << _xc[kj] << DELIMETER;
        }
        fout << endl;
    }

    // Write connectivity: ia_k,0, ia_k,1 etc  in separate lines
    assert( Nelems()*NverticesElements() ==  static_cast<int>(_ia.size()));
    for (int k = 0, kj = 0; k < Nelems(); ++k)
    {
        for (int j = 0; j < NverticesElements(); ++j, ++kj)
        {
            fout << _ia[kj] + OFFSET << DELIMETER;     // C to matlab
        }
        fout << endl;
    }

    // Write vector
    for (int k = 0; k < Nnodes(); ++k)
    {
        fout << v[k] << endl;
    }

    fout.close();
    return;
}


void Mesh::Export_scicomp(std::string const &basename) const
{
    //assert(Nnodes() ==  static_cast<int>(v.size()));  // fits vector length to mesh information?
    string const DELIMETER(" ");    // define the same delimiter as in matlab/ascii_read*.m
    int const    OFFSET(0);
    {
        // Write coordinates into scicomp-file
        string fname(basename + "_coords.txt");
        ofstream fout(fname);                             // open file ASCII mode
        if ( !fout.is_open() )
        {
            cout << "\nFile " << fname << " has not been opened.\n\n" ;
            assert( fout.is_open() && "File not opened."  );
        }

        fout << Nnodes() << endl;
        // Write coordinates: x_k, y_k   in separate lines
        assert( Nnodes()*Ndims() ==  static_cast<int>(_xc.size()));
        for (int k = 0, kj = 0; k < Nnodes(); ++k)
        {
            for (int j = 0; j < Ndims(); ++j, ++kj)
            {
                fout << _xc[kj] << DELIMETER;
            }
            fout << endl;
        }
        fout.close();

    }

    {
        // Write elements into scicomp-file
        string fname(basename + "_elements.txt");
        ofstream fout(fname);                             // open file ASCII mode
        if ( !fout.is_open() )
        {
            cout << "\nFile " << fname << " has not been opened.\n\n" ;
            assert( fout.is_open() && "File not opened."  );
        }

        fout << Nelems() << endl;

        // Write connectivity: ia_k,0, ia_k,1 etc  in separate lines
        assert( Nelems()*NverticesElements() ==  static_cast<int>(_ia.size()));
        for (int k = 0, kj = 0; k < Nelems(); ++k)
        {
            for (int j = 0; j < NverticesElements(); ++j, ++kj)
            {
                fout << _ia[kj] + OFFSET << DELIMETER;     // C to matlab
            }
            fout << endl;
        }
        fout.close();
    }

    return;
}


void Mesh::Visualize_matlab(vector<double> const &v) const
{
    // define external command
    const string exec_m("matlab -nosplash < visualize_results.m");                 // Matlab
    //const string exec_m("octave --no-window-system --no-gui visualize_results.m"); // Octave
    //const string exec_m("flatpak run org.octave.Octave visualize_results.m");      // Octave (flatpak): desktop GH

    const string fname("uv.txt");
    Write_ascii_matlab(fname, v);

    int ierror = system(exec_m.c_str());                                 // call external command

    if (ierror != 0)
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
  
    reNumberEntries(old2new, _bedges);

    reNumberEntries(old2new, _edges);
    sortAscending_2(_edges);       // ascending order of vertices in edge
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

    assert( Nelems()*NverticesElements() ==  static_cast<int>(_ia.size()));
    for (int k = 0, kj = 0; k < Nelems(); ++k)
    {
        fout << 3 << DELIMETER;          // triangular patches
        for (int j = 0; j < NverticesElements()-1; ++j, ++kj)
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



std::vector<int> Mesh::Index_DirichletNodes() const
{
    vector<int> idx(_bedges);                              // copy

    sort(idx.begin(), idx.end());                         // sort
    idx.erase( unique(idx.begin(), idx.end()), idx.end() ); // remove duplicate data

    return idx;
}

// GH
//  only correct for simplices

/*
void Mesh::DeriveEdgeFromVertexBased_fast_2()
{
    assert(NedgesElements() == 4);
    assert(NverticesElements() == 4);   // 3 vertices, 3 edges per element are assumed

    // Store indices of all elements connected to a vertex
    vector<vector<int>> vertex2elems(_nnode, vector<int>(0));
    for (int k = 0; k < Nelems(); ++k)
    {
        for (int i = 0; i < 4; ++i)
        {
            vertex2elems[_ia[4 * k + i]].push_back(k);
        }
    }
    size_t max_neigh = 0;             // maximal number of elements per vertex
    for (auto const &v : vertex2elems)
    {
        max_neigh = max(max_neigh, v.size());
    }
    //cout << endl << vertex2elems << endl;

    // assign edges to elements
    _ea.clear();                      // old data still in _ea without clear()
    _ea.resize(NedgesElements()*Nelems(), -1);
    // Derive the edges
    _edges.clear();
    _nedge = 0;

    // convert also boundary edges
    unsigned int mbc(_bedges.size() / 2);     // number of boundary edges
    _ebedges.clear();
    _ebedges.resize(mbc, -1);
    vector<bool> bdir(_nnode, false);         // vector indicating boundary nodes
    for (size_t kb = 0; kb < _bedges.size(); ++kb )
    {
        bdir.at(_bedges[kb]) = true;
    }

    vector<int> vert_visited;             // already visisted neighboring vertices of k
    vert_visited.reserve(max_neigh);      // avoids multiple (re-)allocations
    for (int k = 0; k < _nnode; ++k)              // vertex k
    {
        vert_visited.clear();
        auto const &elems = vertex2elems[k];      // element neighborhood
        int kedges = static_cast<int>(_edges.size()) / 2; // #edges before vertex k is investigated
        //cout << elems << endl;
// GH: problem, shared edges appear twice.
        int nneigh = elems.size();
        for (int ne = 0; ne < nneigh; ++ne)       // iterate through neighborhood
        {
            int e = elems[ne];                    // neighboring element e
            //cout << "e = " << e << endl;
            for (int i = 4 * e + 0; i < 4 * e + _nvert_e; ++i)   // vertices of element e
            {
                int const vert = _ia[i];
                //cout << "vert: " << vert << "  "<< k << endl;
                if ( vert > k )
                {
                    int ke = -1;
                    auto const iv = find(vert_visited.cbegin(), vert_visited.cend(), vert);
                    if (iv == vert_visited.cend())   // vertex not yet visited
                    {
                        vert_visited.push_back(vert);  // now, vertex vert is visited
                        _edges.push_back(k);           // add the new edge k->vert
                        _edges.push_back(vert);

                        ke = _nedge;
                        ++_nedge;
                        // Is edge ke also a boundary edge?
                        if (bdir[k] && bdir[vert])
                        {
                            size_t kb = 0;
                            while (kb < _bedges.size() && (!( (_bedges[kb] == k && _bedges[kb + 1] == vert) || (_bedges[kb] == vert && _bedges[kb + 1] == k) )) )
                            {
                                kb += 2;
                            }
                            if (kb < _bedges.size())
                            {
                                _ebedges[kb / 2] = ke;
                            }
                        }
                    }
                    else
                    {
                        int offset = iv - vert_visited.cbegin();
                        ke = kedges + offset;
                    }
                    // assign that edge to the edges based connectivity of element e
                    auto ip = find_if(_ea.begin() + 4 * e, _ea.begin() + 4 * (e + 1),
                                      [] (int v) -> bool {return v < 0;} );
                    //cout << ip-_ea.begin()+3*e << "  " << *ip << endl;
                    assert(ip != _ea.cbegin() + 4 * (e + 1)); // data error !
                    *ip = ke;
                }
            }
        }
    }

    assert( Mesh::Check_array_dimensions() );
    return;
}
// HG

// GH
//  only correct for simplices
void Mesh::DeriveEdgeFromVertexBased_fast()
{
    assert(NedgesElements() == 4);
    assert(NverticesElements() == 4);   // 3 vertices, 3 edges per element are assumed

    // Store indices of all elements connected to a vertex
    vector<vector<int>> vertex2elems(_nnode, vector<int>(0));
    for (int k = 0; k < Nelems(); ++k)
    {
        for (int i = 0; i < 4; ++i)
        {
            vertex2elems[_ia[4 * k + i]].push_back(k);
        }
    }
    size_t max_neigh = 0;             // maximal number of elements per vertex
    for (auto const &v : vertex2elems)
    {
        max_neigh = max(max_neigh, v.size());
    }
    //cout << endl << vertex2elems << endl;

    // assign edges to elements
    _ea.clear();                      // old data still in _ea without clear()
    _ea.resize(NedgesElements()*Nelems(), -1);
    // Derive the edges
    _edges.clear();
    _nedge = 0;
    vector<int> vert_visited;             // already visisted neighboring vertices of k
    vert_visited.reserve(max_neigh);      // avoids multiple (re-)allocations
    for (int k = 0; k < _nnode; ++k)              // vertex k
    {
        vert_visited.clear();
        auto const &elems = vertex2elems[k];      // element neighborhood
        int kedges = static_cast<int>(_edges.size()) / 2; // #edges before vertex k is investigated
        //cout << elems << endl;
// GH: problem, shared edges appear twice.
        int nneigh = elems.size();
        for (int ne = 0; ne < nneigh; ++ne)       // iterate through neighborhood
        {
            int e = elems[ne];                    // neighboring element e
            //cout << "e = " << e << endl;
            for (int i = 3 * e + 0; i < 3 * e + _nvert_e; ++i)   // vertices of element e
            {
                int const vert = _ia[i];
                //cout << "vert: " << vert << "  "<< k << endl;
                if ( vert > k )
                {
                    int ke = -1;
                    auto const iv = find(vert_visited.cbegin(), vert_visited.cend(), vert);
                    if (iv == vert_visited.cend())   // vertex not yet visited
                    {
                        vert_visited.push_back(vert);  // now, vertex vert is visited
                        _edges.push_back(k);           // add the new edge k->vert
                        _edges.push_back(vert);

                        ke = _nedge;
                        ++_nedge;
                    }
                    else
                    {
                        int offset = iv - vert_visited.cbegin();
                        ke = kedges + offset;
                    }
                    // assign that edge to the edges based connectivity of element e
                    auto ip = find_if(_ea.begin() + 3 * e, _ea.begin() + 3 * (e + 1),
                                      [] (int v) -> bool {return v < 0;} );
                    //cout << ip-_ea.begin()+3*e << "  " << *ip << endl;
                    assert(ip != _ea.cbegin() + 3 * (e + 1)); // data error !
                    *ip = ke;
                }
            }
        }
    }

    // convert also boundary edges
    unsigned int mbc(_bedges.size() / 2);     // number of boundary edges
    _ebedges.clear();
    _ebedges.resize(mbc, -1);
    for (unsigned int kb = 0; kb < mbc; ++kb)
    {
        int const v1 = min(_bedges[2 * kb], _bedges[2 * kb + 1]); // vertices
        int const v2 = max(_bedges[2 * kb], _bedges[2 * kb + 1]);

        size_t e = 0;
        //   ascending vertex indices for each edge e in _edges
        while (e < _edges.size() && (_edges[e] != v1 || _edges[e + 1] != v2) )
        {
            e += 2;                               // next edge
        }
        assert(e < _edges.size());                // error: no edge found
        _ebedges[kb] = e / 2;                     // index of edge
    }


    assert( Mesh::Check_array_dimensions() );
    return;
}
// HG


#include <utility>             // pair

void Mesh::DeriveEdgeFromVertexBased_slow()
{
    assert(NedgesElements() == 3);
    assert(NverticesElements() == 3);   // 3 vertices, 3 edges per element are assumed

    _ea.resize(NedgesElements()*Nelems());
    vector< pair<int, int> >  edges(0);
    int nedges = 0;

    for (int k = 0; k < Nelems(); ++k)
    {
        array < int, 3 + 1 > ivert{{ _ia[3 * k], _ia[4 * k + 1], _ia[4 * k + 2], _ia[4 * k] }};

        for (int i = 0; i < 3; ++i)
        {
            pair<int, int> e2;          // this edge
            if (ivert[i] < ivert[i + 1])   // guarantee ascending order
            {
                e2 = make_pair(ivert[i], ivert[i + 1]);
            }
            else
            {
                e2 = make_pair(ivert[i + 1], ivert[i]);
            }

            int eki(-1);                // global index of this edge
            auto ip = find(edges.cbegin(), edges.cend(), e2);
            if ( ip == edges.cend() )   // edge not found ==> add that edge
            {
                //cout << "found edge\n";
                edges.push_back(e2);    // add the new edge
                eki = nedges;           // index of this new edge
                ++nedges;

            }
            else
            {
                eki = ip - edges.cbegin(); // index of the edge found
            }
            _ea[3 * k + i] = eki;          // set edge index in edge based connectivity
        }
    }

    assert( nedges == static_cast<int>(edges.size()) );
    _nedge = nedges;                    // set the member variable for number of edges
    _edges.resize(2 * nedges);          // allocate memory for edge storage
    for (int k = 0; k < nedges; ++k)
    {
        _edges[2 * k    ] = edges[k].first;
        _edges[2 * k + 1] = edges[k].second;
    }

    // convert also boundary edges
    unsigned int mbc(_bedges.size() / 2);     // number of boundary edges
    //cout << "AA  " << mbc << endl;
    _ebedges.resize(mbc);
    for (unsigned int kb = 0; kb < mbc; ++kb)
    {
        const auto vv1 = make_pair(_bedges[2 * kb  ], _bedges[2 * kb + 1]); // both
        const auto vv2 = make_pair(_bedges[2 * kb + 1], _bedges[2 * kb  ]); //  directions of edge
        auto ip1 = find(edges.cbegin(), edges.cend(), vv1);
        if (ip1 == edges.cend())
        {
            ip1 = find(edges.cbegin(), edges.cend(), vv2);
            assert(ip1 != edges.cend());          // stop because inconsistency (boundary edge has to be included in edges)
        }
        _ebedges[kb] = ip1 - edges.cbegin();      // index of edge
    }

    assert( Mesh::Check_array_dimensions() );
    return;
}

void Mesh::DeriveVertexFromEdgeBased()
{
    assert(NedgesElements() == 3);
    assert(NverticesElements() == 3);   // 3 vertices, 3 edges per element are assumed

    _ia.resize(NedgesElements()*Nelems()); // NN

    for (int k = 0; k < Nelems(); ++k)
    {
        //vector<int> ivert(6);           // indices of vertices
        array<int, 6> ivert;          // indices of vertices
        for (int j = 0; j < 3; ++j)     // local edges
        {
            int const iedg = _ea[3 * k + j]; // index of one edge in triangle
            ivert[2 * j  ] = _edges[2 * iedg  ]; // first  vertex of edge
            ivert[2 * j + 1] = _edges[2 * iedg + 1]; // second vertex of edge
        }
        sort(ivert.begin(), ivert.end()); // unique indices are needed
        auto const ip = unique(ivert.begin(), ivert.end());
        assert( ip - ivert.begin() == 3 );
        for (int i = 0; i < 3; ++i)     // vertex based element connectivity
        {
            _ia[3 * k + i] = ivert[i];
        }
    }

    // convert also boundary edges
    unsigned int mbc(_ebedges.size());       // number of boundary edges
    _bedges.resize(2 * mbc);
    for (unsigned int k = 0; k < mbc; ++k)
    {
        const auto ke = _ebedges[k];        // edge index
        _bedges[2 * k  ] = _edges[2 * ke  ];
        _bedges[2 * k + 1] = _edges[2 * ke + 1];
    }


    return;
}
*/














vector<vector<int>> Mesh::Node2NodeGraph_2() const
{
    vector<vector<int>> v2v(_nnode, vector<int>(0));    // stores the vertex to vertex connections

    ////--------------
    vector<int> cnt(_nnode);
    for (size_t i = 0; i < _ia.size(); ++i)  ++cnt[_ia[i]]; // determine number of entries per vertex
    for (size_t k = 0; k < v2v.size(); ++k)
    {
        v2v[k].resize(_nvert_e * cnt[k]);              //    and allocate the memory for that vertex
        cnt[k] = 0;
    }
    ////--------------

    for (int e = 0; e < _nelem; ++e)
    {
        int const basis = e * _nvert_e;                  // start of vertex connectivity of element e
        for (int k = 0; k < _nvert_e; ++k)
        {
            int const v = _ia[basis + k];
            for (int l = 0; l < _nvert_e; ++l)
            {
                v2v[v][cnt[v]] = _ia[basis + l];
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

vector<vector<int>> Mesh::Node2NodeGraph_1() const
{
    vector<vector<int>> v2v(_nnode, vector<int>(0));    // stores the vertex to vertex connections

    for (int e = 0; e < _nelem; ++e)
    {
        int const basis = e * _nvert_e;                  // start of vertex connectivity of element e
        for (int k = 0; k < _nvert_e; ++k)
        {
            int const v = _ia[basis + k];
            for (int l = 0; l < _nvert_e; ++l)
            {
                v2v[v].push_back(_ia[basis + l]);
            }
        }
    }
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


Mesh::Mesh(std::string const &fname)
    : Mesh(3, 4, 4) // two dimensions, 3 vertices, 3 dofs, 3 edges per element
{
    ReadVertexBasedMesh(fname);
    //DeriveEdgeFromVertexBased();        // Generate also the edge based information
    
     //// GH: Check permuted numbering
    //vector<int> perm(Nnodes());
    //iota(rbegin(perm),rend(perm),0);
    //random_shuffle(begin(perm),end(perm));
    //PermuteVertices(perm);   
    //cout << " P E R M U T E D !" << endl;
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
    assert(ndim == 3 && nvert_e == 4);

    // Allocate memory
    Resize_Coords(nnode, ndim);                 // coordinates in 2D [nnode][ndim]
    Resize_Connectivity(nelem, nvert_e);        // connectivity matrix [nelem][nvert]

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

    // additional read of boundary information (only start/end point)
    int nbedges;
    ifs >> nbedges;

    _bedges.resize(nbedges * 2);
    for (int k = 0; k < nbedges * 2; ++k)
    {
        ifs >> _bedges[k];
        _bedges[k] -= OFFSET;            // Matlab to C indexing
    }

    return;
}


double Mesh::largestAngle(int idx) const
{
    assert( 0<=idx && idx<Nelems() );
    assert( 2==Ndims() );
    int const NE=3;
    
    array<int,3> const gvert{_ia[NE*idx],_ia[NE*idx+1],_ia[NE*idx+2]};
    array<double,NE> vcos;
    for (int vm=0; vm<NE; ++vm)
    {
        int const vl = (vm-1+NE) % NE;
        int const vr = (vm+1   ) % NE;
        array<double,2> vec_l{ _xc[2*gvert[vl]]-_xc[2*gvert[vm]], _xc[2*gvert[vl]+1]-_xc[2*gvert[vm]+1] };
        array<double,2> vec_r{ _xc[2*gvert[vr]]-_xc[2*gvert[vm]], _xc[2*gvert[vr]+1]-_xc[2*gvert[vm]+1] };
        vcos[vm] = (vec_l[0]*vec_r[0]+vec_l[1]*vec_r[1])/hypot(vec_l[0],vec_l[1])/hypot(vec_r[0],vec_r[1]);
    }
    double vmax = *min_element(cbegin(vcos),cend(vcos));
    
    return acos(vmax);
}

vector<double> Mesh::getLargestAngles() const
{
    vector<double> angles(Nelems());
    for (int elem=0; elem<Nelems(); ++elem)
    {
        angles[elem] = largestAngle(elem);
    }
    return angles;
}


bool Mesh::checkObtuseAngles() const
{
    vector<double> const angles=getLargestAngles();    
    // C++20: #include <numbers>  std::numbers::pi or std::numbers::pi_v<double>
    bool bb=any_of(cbegin(angles),cend(angles), [](double a) -> bool {return a>M_PI/2.0;} );
    
    if (bb)     // find elements with largest angles
    {   
        cout << "!!  Those elements with largest inner angle > pi/2  !!\n";
        //vector<double> lpi2(size(angles));
        //auto ip=copy_if(cbegin(angles),end(angles),begin(lpi2), [](double a) -> bool {return a>M_PI/2.0;} );
        //lpi2.erase(ip,end(lpi2));
        //cout << "Lang: " << size(lpi2) << endl;
        vector<int> large = sort_indexes_desc(angles);
        size_t num_show = min(10ul,size(large));
        for (size_t k=0; k<num_show; ++k)
        {
            if (angles[large[k]] > M_PI/2.0)
            {
                cout << "elem [" << large[k] << "] : " << angles[large[k]] << endl;
            }
        }
    }
    
    return bb;
}











bool Mesh::Check_array_dimensions() const
{
    bool b_ia = static_cast<int>(_ia.size() / _nvert_e) == _nelem;
    if (!b_ia)  cerr << "misfit: _nelem vs. _ia" << endl;

    bool b_xc = static_cast<int>(_xc.size() / _ndim) == _nnode;
    if (!b_xc)  cerr << "misfit: _nnode vs. _xc" << endl;

    bool b_ea = static_cast<int>(_ea.size() / _nedge_e) == _nelem;
    if (!b_ea)  cerr << "misfit: _nelem vs. _ea" << endl;

    bool b_ed = static_cast<int>(_edges.size() / 2) == _nedge;
    if (!b_ed)  cerr << "misfit: _nedge vs. _edges" << endl;


    return b_ia && b_xc && b_ea && b_ed;
}

void Mesh::Del_EdgeConnectivity()
{
    _nedge = 0;              //!< number of edges in mesh
    _edges.resize(0);        //!< edges of mesh (vertices ordered ascending)
    _edges.shrink_to_fit();
    _ea.resize(0);           //!< edge based element connectivity
    _ea.shrink_to_fit();
    _ebedges.resize(0);      //!< boundary edges [nbedges]
    _ebedges.shrink_to_fit();
    return;
}



// ####################################################################

RefinedMesh::RefinedMesh(Mesh const &cmesh, std::vector<bool> const &ibref)
//: Mesh(cmesh), _cmesh(cmesh), _ibref(ibref), _nref(0), _vfathers(0)
    : Mesh(cmesh), _ibref(ibref), _nref(0), _vfathers(0)
{
    if (_ibref.size() == 0)                  // refine all elements
    {
        //
        RefineAllElements();
    }
    else
    {
        cout << endl << "  Adaptive Refinement not implemented yet." << endl;
        assert(_ibref.size() != 0);
    }
}

RefinedMesh::~RefinedMesh()
{}

Mesh RefinedMesh::RefineElements(std::vector<bool> const & /*ibref*/)
{
    Mesh new_mesh(_ndim, _nvert_e, _ndof_e);
    cout << " NOT IMPLEMENTED: Mesh::RefineElements" << endl;

////  initialize new coorsinates with the old one
    //auto new_coords = new_mesh.GetCoords();
    //new_coords = _xc;                      // copy coordinates from old mesh

//// access vertex connectivite, edge connectiviy and edge information of new mesh
    //auto new_ia = new_mesh.GetConnectivity();
    //auto new_ea = new_mesh.GetEdgeConnectivity();
    //auto new_edges = new_mesh.GetEdges();

////  storing the parents of edges and vertices


    //assert( new_ia.size()== new_ea.size() );
    //new_mesh.SetNnode( new_coords.size() );
    //new_mesh.SetNelem( new_ia.size()/3 );
    //new_mesh._nedge = new_edges.size()/2;

    return new_mesh;
}

//JF


/*
void RefinedMesh::RefineAllElements(int nref)
{
    cout << "\n############   Refine Mesh " << nref << " times ";
    auto tstart = clock();
    //DeriveEdgeFromVertexBased();          // ensure that edge information is available

    for (int kr = 0; kr < nref; ++kr)
    {
        //DeriveEdgeFromVertexBased();          // ensure that edge information is available // GH: not needed in each loop

        auto old_ea(_ea);                     // save old edge connectivity
        auto old_edges(_edges);               // save old edges
        auto old_nedges(Nedges());
        auto old_nnodes(Nnodes());
        auto old_nelems(Nelems());

        //  the new vertices will be appended to the coordinates in _xc

        vector<int> edge_sons(2 * old_nedges); // 2 sons for each edge

        //   --  Derive the fine edges ---
        int new_nedge = 2 * old_nedges + 3 * old_nelems; // #edges in new mesh
        int new_nelem = 4 * old_nelems;     // #elements in new mesh
        int new_nnode = old_nnodes + old_nedges; // #nodes in new mesh

        _xc.reserve(2 * new_nnode);
        // store the 2 fathers of each vertex (equal fathers denote original coarse vertex)
        _vfathers.resize(2 * old_nnodes);
        for (int vc = 0; vc < old_nnodes; ++vc)
        {
            _vfathers[2 * vc  ] = vc;    // equal fathers denote original coarse vertex
            _vfathers[2 * vc + 1] = vc;
        }

        _ia.clear();
        _ea.clear();
        _ea.resize(new_nelem * 3);
        _edges.clear();
        _edges.resize(2 * new_nedge);          // vertices of edges [v_0, v_1;v_0, v_1; ...]
        vector<int> e_son(2 * old_nedges); // sons of coarse edges [s_0, s_1; s_0, s_1; ...]

        // split all coarse edges and append the new nodes
        int kf = 0;                            // index of edges in fine mesh
        int vf = old_nnodes;              // index of new vertex in fine grid
        for (int kc = 0; kc < old_nedges; ++kc)   // index of edges in coarse mesh
        {
            //
            int v1 = old_edges[2 * kc];        // vertices of old edge
            int v2 = old_edges[2 * kc + 1];
            // append coordinates of new vertex
            double xf = 0.5 * ( _xc[2 * v1  ] + _xc[2 * v2  ] );
            double yf = 0.5 * ( _xc[2 * v1 + 1] + _xc[2 * v2 + 1] );
            _xc.push_back(xf);
            _xc.push_back(yf);
            // fathers of vertex  vf
            _vfathers.push_back(v1);
            _vfathers.push_back(v2);

            // split old edge into two edges
            _edges[2 * kf    ] = v1;             // coarse vertex 1
            _edges[2 * kf + 1] = vf;             //   to new fine vertex
            e_son[2 * kc    ]  = kf;             // son edge
            ++kf;
            _edges[2 * kf    ] = vf;             // new fine vertex
            _edges[2 * kf + 1] = v2;             //   to coarse vertex 2
            e_son[2 * kc + 1]  = kf;             // son edge

            ++vf;
            ++kf;
        }
        _xc.shrink_to_fit();
        _vfathers.shrink_to_fit();

        // -- derive the fine mesh elements --
        //    creates additional fine edges

        for (int kc = 0; kc < old_nelems; ++kc)   // index of elements in coarse mesh
        {
            array<array<int, 3>, 3 * 2> boundary; // fine scale vertices and edges as boundary of old element
            //boundary[ ][0], boundary[ ][1] ..vertices boundary[ ][2] edge

            for (int j = 0; j < 3; ++j)           // each edge in element
            {
                int ce = old_ea[3 * kc + j];      // coarse edge number

                int s1 = e_son[2 * ce    ];       // son edges of that coarse edge
                int s2 = e_son[2 * ce + 1];
                boundary[2 * j][2] = s1;          //add boundary edge
                boundary[2 * j][0] = _edges[2 * s1 + 0];
                boundary[2 * j][1] = _edges[2 * s1 + 1];
                if (boundary[2 * j][0] > boundary[2 * j][1]) swap(boundary[2 * j][0], boundary[2 * j][1]); 		// fine vertices always in 2nd entry
                boundary[2 * j + 1][2] = s2;      //add boundary edge
                boundary[2 * j + 1][0] = _edges[2 * s2 + 0];
                boundary[2 * j + 1][1] = _edges[2 * s2 + 1];
                if (boundary[2 * j + 1][0] > boundary[2 * j + 1][1]) swap(boundary[2 * j + 1][0], boundary[2 * j + 1][1]);
            }

            sort(boundary.begin(), boundary.end());		// sort -> edges with same coarse vertex will be neighbors

            int interior_1 = 2 * old_nedges + kc * 3;	// add interior edges
            int interior_2 = 2 * old_nedges + kc * 3 + 1;
            int interior_3 = 2 * old_nedges + kc * 3 + 2;

            _edges[interior_1 * 2    ] = boundary[0][1]; // add interior edges
            _edges[interior_1 * 2 + 1] = boundary[1][1];

            _edges[interior_2 * 2    ] = boundary[2][1];
            _edges[interior_2 * 2 + 1] = boundary[3][1];

            _edges[interior_3 * 2    ] = boundary[4][1];
            _edges[interior_3 * 2 + 1] = boundary[5][1];

            _ea[kc * 3 * 4    ] = boundary[0][2];       // add 4 new elements with 3 edges for every old element
            _ea[kc * 3 * 4 + 1] = boundary[1][2];
            _ea[kc * 3 * 4 + 2] = interior_1;

            _ea[kc * 3 * 4 + 3] = boundary[2][2];
            _ea[kc * 3 * 4 + 4] = boundary[3][2];
            _ea[kc * 3 * 4 + 5] = interior_2;

            _ea[kc * 3 * 4 + 6] = boundary[4][2];
            _ea[kc * 3 * 4 + 7] = boundary[5][2];
            _ea[kc * 3 * 4 + 8] = interior_3;

            _ea[kc * 3 * 4 + 9] = interior_1;
            _ea[kc * 3 * 4 + 10] = interior_2;
            _ea[kc * 3 * 4 + 11] = interior_3;
        }

// GH: ToDo:  _bedges  has to updated for the new mesh //!< boundary edges [nbedges][2] storing start/end vertex
//     Pass the refinement information to the boundary edges (edge based)
        auto old_ebedges(_ebedges);           // save original boundary edges [nbedges] (edge based storage)
        unsigned int old_nbedges(old_ebedges.size());

        _ebedges.resize(2 * old_nbedges);    // each old boundary edge will be bisected
        unsigned int kn = 0;                 // index of new boundary edges
        for (unsigned int ke = 0; ke < old_nbedges; ++ke)   // index of old boundary edges
        {
            const auto kc = old_ebedges[ke];
            _ebedges[kn] = e_son[2 * kc    ];
            ++kn;
            _ebedges[kn] = e_son[2 * kc + 1];
            ++kn;
        }
// HG
        // set new mesh parameters
        SetNelem(new_nelem);
        SetNnode(new_nnode);
        SetNedge(new_nedge);

        {
// Cuthill-McKee reordering
//     Increases mesh generation time by factor 5 -  but solver is faster.
            auto const perm = cuthill_mckee_reordering(_edges);
            PermuteVertices_EdgeBased(perm);
        }

        //DeriveVertexFromEdgeBased();
        assert( RefinedMesh::Check_array_dimensions() );

        ++_nref;                            // track the number of refinements
    }

    double duration = static_cast<double>(clock() - tstart) / CLOCKS_PER_SEC;
    cout << "finished in  " <<  duration  << " sec.    ########\n";

    return;
}


void Mesh::PermuteVertices_EdgeBased(vector<int> const &old2new)
{
//      permute vertices _edges
    auto const edges_old(_edges);
    for (size_t k = 0; k < _edges.size(); k += 2)
    {
        _edges[k    ] = old2new[edges_old[k    ]];
        _edges[k + 1] = old2new[edges_old[k + 1]];
        if (_edges[k] > _edges[k + 1])
            swap(_edges[k], _edges[k + 1]);
    }
//      permute coordinates
    auto const coord_old(_xc);
    for (size_t k = 0; k < _xc.size() / 2; ++k)
    {
        _xc[2 * old2new[k]    ] = coord_old[2 * k    ];
        _xc[2 * old2new[k] + 1] = coord_old[2 * k + 1];
    }
    return;
}
*/

void RefinedMesh::PermuteVertices_EdgeBased(vector<int> const &old2new)
{
    Mesh::PermuteVertices_EdgeBased(old2new);
//      permute fathers of a vertex
    auto const old_fathers(_vfathers);
    for (size_t k = 0; k < _vfathers.size() / 2; ++k)
    {
        _vfathers[2 * old2new[k]    ] = old_fathers[2 * k    ];
        _vfathers[2 * old2new[k] + 1] = old_fathers[2 * k + 1];
    }
    return;
}


bool RefinedMesh::Check_array_dimensions() const
{
    const bool bp = Mesh::Check_array_dimensions();

    const bool bvf = (static_cast<int>(_vfathers.size()) / 2 == Nnodes());

    return bp && bvf;

}
// #####################################################################


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




/*






// #####################################################################
Mesh_2d_3_square::Mesh_2d_3_square(int nx, int ny, int myid, int procx, int procy)
    : Mesh(2, 3, 3, 3), // two dimensions, 3 vertices, 3 dofs, 3 edges per element
      _myid(myid), _procx(procx), _procy(procy), _neigh{{ -1, -1, -1, -1}}, _color(0),
_xl(0.0), _xr(1.0), _yb(0.0), _yt(1.0), _nx(nx), _ny(ny)
{
    //void IniGeom(int const myid, int const procx, int const procy, int neigh[], int &color)
    int const ky = _myid / _procx;
    int const kx = _myid % _procy;	//    MOD(myid,procx)
    // Determine the neighbors of domain/rank myid
    _neigh[0] = (ky == 0)       ?  -1 : _myid - _procx;    //   South
    _neigh[1] = (kx == _procx - 1) ?  -1 : _myid + 1;      //   East
    _neigh[2] = (ky == _procy - 1) ?  -1 : _myid + _procx; //   North
    _neigh[3] = (kx == 0)       ?  -1 : _myid - 1;         //   West

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

    return;
}

Mesh_2d_3_square::~Mesh_2d_3_square()
{}


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
    int const start[4] = { bl, br, tl, bl};
    int const   end[4] = { br, tr, tr, tl};
    int const  step[4] = { dx, dy, dx, dy};

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
    sort(idx.begin(), idx.end());                          // sort
    idx.erase( unique(idx.begin(), idx.end()), idx.end() ); // remove duplicate data

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
}*/

// ####################################################################

