%% Visualize results
%
%   flatpak run org.octave.Octave <filename>
%      or
%   octave --no-window-system --no-gui  -qf <filename>
%
%      or
%   matlab -nosplash <   <filename>

clear all
clc

%%
fname = 'uv.txt';

[xc,ia,v] = ascii_read_meshvector(fname);

h = trisurf(ia, xc(:,1), xc(:,2), v);

waitfor(h)                     % wait for closing the figure