% Square: 
%   flatpak run org.octave.Octave <filename>
%      or
%   octave --no-window-system --no-gui  -qf <filename>

clear all
clc
% %% L-shape
% g=[2 0 2 0 0 1 0;        % #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
%    2 2 2 0 1 1 0;
%    2 2 1 1 0.5 1 0;
%    2 1 1 0.5 2 1 0;
%    2 1 0 2 2 1 0;
%    2 0 0 2 0 1 0]';

%% square
g=[2 0 1 0 0 1 0;        % #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
   2 1 1 0 1 1 0;
   2 1 0 1 1 1 0;
   2 0 0 1 0 1 0]';

%[p,e,t] = initmesh(g,'hmax',0.01); 
[p,e,t] = initmesh(g,'hmax',0.5); 
pdemesh(p,e,t)

%% GH
% output from <https://de.mathworks.com/help/pde/ug/initmesh.html initmesh>
%
% coordinates  p: [2][nnode]
% connectivity t: [4][nelem]   with  t(4,:) are the subdomain numbers
% edges        e: [7][nedges]  boundary edges
%                              e([1,2],:) - start/end vertex of edge
%                              e([3,4],:) - start/end values
%                              e(5,:)     - segment number
%                              e([6,7],:) - left/right subdomain

ascii_write_mesh( p, t, e, mfilename);



% tmp=t(1:3,:)

