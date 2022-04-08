function ascii_write_mesh( xc, ia, e, basename)
%
% Saves the 2D triangular mesh in the minimal way (only coordinates, vertex connectivity, minimal boundary edge info)
%  in an ASCII file.
%  Matlab indexing is stored  (starts with 1).
% 
% The output file format is compatible with Mesh_2d_3_matlab:Mesh_2d_3_matlab(std::string const &fname) in jacobi_oo_stl/geom.h
%
% IN:
% coordinates  xc: [2][nnode]
% connectivity ia: [4][nelem]   with  t(4,:) are the subdomain numbers
% edges         e: [7][nedges]  boundary edges
%                              e([1,2],:) - start/end vertex of edge
%                              e([3,4],:) - start/end values
%                              e(5,:)     - segment number
%                              e([6,7],:) - left/right subdomain
%        basename: file name without extension
% 
% Data have been generated via <https://de.mathworks.com/help/pde/ug/initmesh.html initmesh>.
%
fname = [basename, '.txt'];

nnode = int32(size(xc,2));
ndim  = int32(size(xc,1));
nelem = int32(size(ia,2));
nvert_e = int32(3);


dlmwrite(fname,nnode,'delimiter','\t','precision',16)                 % number of nodes
dlmwrite(fname,ndim,'-append','delimiter','\t','precision',16)        % space dimension
dlmwrite(fname,nelem,'-append','delimiter','\t','precision',16)       % number of elements
dlmwrite(fname,nvert_e,'-append','delimiter','\t','precision',16)     % number of vertices per element

% dlmwrite(fname,xc(:),'-append','delimiter','\t','precision',16)       % coordinates
dlmwrite(fname,xc([1,2],:).','-append','delimiter','\t','precision',16) % coordinates

% no subdomain info transferred
tmp=int32(ia(1:3,:));
% dlmwrite(fname,tmp(:),'-append','delimiter','\t','precision',16)      % connectivity in Matlab indexing
dlmwrite(fname,tmp(:,:).','-append','delimiter','\t','precision',16)    % connectivity in Matlab indexing

% store only start and end point of boundary edges,
nbedges = size(e,2);
dlmwrite(fname,nbedges,'-append','delimiter','\t','precision',16)     % number boundary edges
tmp=int32(e(1:2,:));
% dlmwrite(fname,tmp(:),'-append','delimiter','\t','precision',16)    % boundary edges in Matlab indexing
dlmwrite(fname,tmp(:,:).','-append','delimiter','\t','precision',16)  % boundary edges in Matlab indexing

end
