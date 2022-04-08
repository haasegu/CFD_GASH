function [ xc, ia, v ] = ascii_read_meshvector( fname )
%
% Loads the 2D triangular mesh (coordinates, vertex connectivity) 
%  together with values on its vertices from an ASCII file.
%   Matlab indexing is stored  (starts with 1).
% 
% The input file format is compatible 
%  with Mesh_2d_3_matlab:Write_ascii_matlab(..) in jacobi_oo_stl/geom.h
%
%
%  IN: fname - filename
% OUT: xc    - coordinates
%      ia    - mesh connectivity
%      v     - solution vector

DELIMETER = ' ';

fprintf('Read file  %s\n',fname)

% Read mesh constants
nn = dlmread(fname,DELIMETER,[0 0 0 3]);  %% row_1, col_1, row_2, col_2  in C indexing!!!
nnode = nn(1);
ndim  = nn(2);
nelem = nn(3);
nvert = nn(4);

% Read coordinates
row_start = 0+1;
row_end   = 0+nnode;
xc = dlmread(fname,DELIMETER,[row_start 0 row_end ndim-1]);

% Read connectivity
row_start = row_end+1;
row_end   = row_end+nelem;
ia = dlmread(fname,DELIMETER,[row_start 0 row_end nvert-1]);

% Read solution
row_start = row_end+1;
row_end   = row_end+nnode;
v = dlmread(fname,DELIMETER,[row_start 0 row_end 0]);
end


