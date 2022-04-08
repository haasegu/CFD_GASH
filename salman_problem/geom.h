#ifndef GEOM_FILE
#define GEOM_FILE
#include <array>
#include <functional>             // function; C++11
#include <string>
#include <vector>

/**
 * Basis class for finite element meshes.
 */
class Mesh
{
public:
    /**
     * Constructor initializing the members with default values.
      *
      * @param[in] ndim  space dimensions (dimension for coordinates)
      * @param[in] nvert_e  number of vertices per element (dimension for connectivity)
      * @param[in] ndof_e   degrees of freedom per element (= @p nvert_e for linear elements)
     */
    explicit Mesh(int ndim, int nvert_e = 0, int ndof_e = 0);

    /**
     * Destructor.
     *
     * See clang warning on
     * <a href="https://stackoverflow.com/questions/28786473/clang-no-out-of-line-virtual-method-definitions-pure-abstract-c-class/40550578">weak-vtables</a>.
     */
    virtual ~Mesh();

    /**
     * Number of finite elements in (sub)domain.
     * @return number of elements.
     */
    int Nelems() const
    {
        return _nelem;
    }

    /**
     * Global number of vertices for each finite element.
     * @return number of vertices per element.
     */
    int NverticesElements() const
    {
        return _nvert_e;
    }

    /**
     * Global number of degrees of freedom (dof) for each finite element.
     * @return degrees of freedom per element.
     */
    int NdofsElement() const
    {
        return _ndof_e;
    }

    /**
     * Number of vertices in mesh.
     * @return number of vertices.
     */
    int Nnodes() const
    {
        return _nnode;
    }

    /**
     * Space dimension.
     * @return number of dimensions.
     */
    int Ndims() const
    {
        return _ndim;
    }

    /**
     * (Re-)Allocates memory for the element connectivity and redefines the appropriate dimensions.
     *
     * @param[in] nelem    number of elements
     * @param[in] nvert_e  number of vertices per element
     */
    void Resize_Connectivity(int nelem, int nvert_e)
    {
        SetNelem(nelem);               // number of elements
        SetNverticesElement(nvert_e);  // vertices per element
        _ia.resize(nelem * nvert_e);
    }

    /**
     * Read connectivity information (g1,g2,g3)_i.
     * @return connectivity vector [nelems*ndofs].
     */
    const std::vector<int>  &GetConnectivity() const
    {
        return _ia;
    }

    /**
     * Access/Change connectivity information (g1,g2,g3)_i.
     * @return connectivity vector [nelems*ndofs].
     */
    std::vector<int>  &GetConnectivity()
    {
        return _ia;
    }

    /**
     * (Re-)Allocates memory for the element connectivity and redefines the appropriate dimensions.
     *
     * @param[in] nnodes    number of nodes
     * @param[in] ndim      space dimension
     */
    void Resize_Coords(int nnodes, int ndim)
    {
        SetNnode(nnodes);       // number of nodes
        SetNdim(ndim);          // space dimension
        _xc.resize(nnodes * ndim);
    }

    /**
     * Read coordinates of vertices (x,y)_i.
     * @return coordinates vector [nnodes*2].
     */
    const std::vector<double> &GetCoords() const
    {
        return _xc;
    }

    /**
     * Access/Change coordinates of vertices (x,y)_i.
     * @return coordinates vector [nnodes*2].
     */
    std::vector<double> &GetCoords()
    {
        return _xc;
    }

    /**
     * Calculate values in vector @p v via function @p func(x,y)
     * @param[in] v     vector
     * @param[in] func  function of (x,y) returning a double value.
     */
    void SetValues(std::vector<double> &v, const std::function<double(double, double)> &func) const;

    /**
     * Prints the information for a finite element mesh
     */
    void Debug() const;

    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions
     * @return index vector.
     */
    virtual std::vector<int> Index_DirichletNodes() const = 0;

    /**
     * Write vector @p v toghether with its mesh information to an ASCii file @p fname.
      *
      * The data are written in C-style.
      *
      * @param[in] fname  file name
     * @param[in] v      vector
     */
    void Write_ascii_matlab(std::string const &fname, std::vector<double> const &v) const;

    /**
     * Visualize @p v together with its mesh information via matlab or octave.
      *
      * Comment/uncomment those code lines in method Mesh:Visualize (geom.cpp)
      * that are supported on your system.
      *
      * @param[in] v      vector
      *
      * @warning matlab files ascii_read_meshvector.m  visualize_results.m
      *          must be in the executing directory.
     */
    void Visualize(std::vector<double> const &v) const;
    
    /**
      * All data containing vertex numbering are renumbered or sorted according to
      * the permutation @p permut_old2new .
      *
      * @param[in] old2new      permutation of vertex indices: old2new[k] stores the new index of old index k
      *
     */
    virtual void PermuteVertices(std::vector<int> const& old2new);


protected:
    void SetNelem(int nelem)
    {
        _nelem = nelem;
    }

    void SetNverticesElement(int nvert)
    {
        _nvert_e = nvert;
    }

    void SetNdofsElement(int ndof)
    {
        _ndof_e = ndof;
    }

    void SetNnode(int nnode)
    {
        _nnode = nnode;
    }

    void SetNdim(int ndim)
    {
        _ndim = ndim;
    }

private:
    int _nelem;         //!< number elements
    int _nvert_e;       //!< number of vertices per element
    int _ndof_e;        //!< degrees of freedom (d.o.f.) per element
    int _nnode;         //!< number nodes/vertices
    int _ndim;          //!< space dimension of the problem (1, 2, or 3)
    std::vector<int> _ia;    //!< element connectivity
    std::vector<double> _xc; //!< coordinates
};

/**
 * 2D finite element mesh of the square consiting of linear triangular elements.
 */
class Mesh_2d_3_square: public Mesh
{
public:
    /**
     * Generates the f.e. mesh for the unit square.
     *
     * @param[in] nx    number of discretization intervals in x-direction
     * @param[in] ny    number of discretization intervals in y-direction
     * @param[in] myid  my MPI-rank / subdomain
     * @param[in] procx number of ranks/subdomains in x-direction
     * @param[in] procy number of processes in y-direction
    */
    Mesh_2d_3_square(int nx, int ny, int myid = 0, int procx = 1, int procy = 1);

    /**
     * Destructor
     */
    ~Mesh_2d_3_square() override
    {}

    /**
     * Set solution vector based on a tensor product grid in the rectangle.
     * @param[in] u solution vector
     */
    void SetU(std::vector<double> &u) const;

    /**
     * Set right hand side (rhs) vector on a tensor product grid in the rectangle.
     * @param[in] f rhs vector
     */
    void SetF(std::vector<double> &f) const;

    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions
     * @return index vector.
     */
    std::vector<int> Index_DirichletNodes() const override;

    /**
      * Stores the values of vector @p u of (sub)domain into a file @p name for further processing in gnuplot.
      * The file stores rowise the x- and y- coordinates together with the value from  @p u .
      * The domain [@p xl, @p xr] x [@p yb, @p yt] is discretized into @p nx x @p ny intervals.
      *
      * @param[in] name  basename of file name (file name will be extended by the rank number)
      * @param[in] u     local vector
      *
      * @warning   Assumes tensor product grid in unit square; rowise numbered
      *            (as generated in class constructor).
      *            The output is provided for tensor product grid visualization
      *            ( similar to Matlab-surf() ).
      *
      * @see Mesh_2d_3_square
      */
    void SaveVectorP(std::string const &name, std::vector<double> const &u) const;

    // here will still need to implement in the class
    //  GetBound(), AddBound()
    //  or better a generalized way with indices and their appropriate ranks for MPI communication

private:
    /**
      * Determines the coordinates of the dicretization nodes of the domain [@p xl, @p xr] x [@p yb, @p yt]
      * which is discretized into @p nx x @p ny intervals.
      * 
      * @param[in] nx   number of discretization intervals in x-direction
      * @param[in] ny    number of discretization intervals in y-direction
      * @param[in] xl    x-coordinate of left boundary
      * @param[in] xr    x-coordinate of right boundary
      * @param[in] yb    y-coordinate of lower boundary
      * @param[in] yt    y-coordinate of upper boundary
      * @param[out] xc   coordinate vector of length 2n with x(2*k,2*k+1) as coodinates of node k
      */

    void GetCoordsInRectangle(int nx, int ny, double xl, double xr, double yb, double yt,
                              double xc[]);
    /**
      * Determines the element connectivity of linear triangular elements of a FEM discretization
      * of a rectangle using @p nx x @p ny equidistant intervals for discretization.
      * @param[in] nx    number of discretization intervals in x-direction
      * @param[in] ny    number of discretization intervals in y-direction
      * @param[out] ia   element connectivity matrix with ia(3*s,3*s+1,3*s+2) as node numbers od element s
      */
    void GetConnectivityInRectangle(int nx, int ny, int ia[]);

private:
    int _myid;          //!< my MPI rank
    int _procx;         //!< number of MPI ranks in x-direction
    int _procy;         //!< number of MPI ranks in y-direction
    std::array<int, 4> _neigh; //!< MPI ranks of neighbors (negative: no neighbor but b.c.)
    int _color;         //!< red/black coloring (checker board) of subdomains

    double _xl;         //!< x coordinate of lower left  corner of square
    double _xr;         //!< x coordinate of lower right corner of square
    double _yb;         //!< y coordinate or lower left  corner of square
    double _yt;         //!< y coordinate of upper right corner of square
    int    _nx;         //!< number of intervals in x-direction
    int    _ny;         //!< number of intervals in y-direction
};

// #################### still some old code (--> MPI) ############################
/**
 * Copies the values of @p w corresponding to boundary @p ib
 * onto vector s.  South (ib==1), East (ib==2), North (ib==3), West (ib==4).
 * The vector @p s has to be long enough!!
 * @param[in] ib    my local boundary
 * @param[in] nx    number of discretization intervals in x-direction
 * @param[in] ny    number of discretization intervals in y-direction
 * @param[in] w     vector for all nodes of local discretization
 * @param[out] s    short vector with values on boundary @p ib
*/
// GH_NOTE: Absicherung bei s !!
void GetBound(int ib, int nx, int ny, double const w[], double s[]);


/**
 * Computes @p w := @p w + @p s  at the interface/boundary nodes on the
 * boundary @p ib .  South (ib==1), East (ib==2), North (ib==3), West (ib==4)
 * @param[in] ib    my local boundary
 * @param[in] nx    number of discretization intervals in x-direction
 * @param[in] ny    number of discretization intervals in y-direction
 * @param[in,out] w vector for all nodes of local discretization
 * @param[in] s     short vector with values on boundary @p ib
*/
void AddBound(int ib, int nx, int ny, double w[], double const s[]);

// #################### Mesh from Matlab ############################
/**
 * 2D finite element mesh of the square consiting of linear triangular elements.
 */
class Mesh_2d_3_matlab: public Mesh
{
public:
    /**
     * Reads mesh data from a binary file.
     *
     * File format, see ascii_write_mesh.m
     *
     * @param[in] fname file name
    */
    explicit Mesh_2d_3_matlab(std::string const &fname);

    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions.
     * @return index vector.
      *
      * @warning All boundary nodes are considered as Dirchlet nodes.
     */
    std::vector<int> Index_DirichletNodes() const override;
    
    /**
      * All data containing vertex numbering are renumbered or sorted according to
      * the permutation @p permut_old2new .
      *
      * @param[in] old2new      permutation of vertex indices: old2new[k] stores the new index of old index k
      *
     */
    void PermuteVertices(std::vector<int> const& old2new) override;

private:
    /**
     * Determines the indices of those vertices with Dirichlet boundary conditions
     * @return index vector.
     */
    int Nnbedges() const
    {
        return static_cast<int>(bedges.size());
    }

    std::vector<int> bedges;     //!< boundary edges [nbedges][2] storing start/end vertex

};

#endif
