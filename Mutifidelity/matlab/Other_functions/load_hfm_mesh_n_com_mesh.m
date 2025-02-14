function [Mesh_com,Mesh_high,h_high,QoI_h] = load_hfm_mesh_n_com_mesh(com_mesh_level,hfm_mesh_level)


% common mesh
n = 10*2^com_mesh_level;
N = 2*n+1; % #of mesh grids
Mesh_com = linspace(0,2,N)';
% h = 2/(2*n);
% QoI_l = N;


%high-didelity
n = 10*2^hfm_mesh_level;
N = 2*n+1; % #of mesh grids
Mesh_high = linspace(0,2,N)';
h_high = 2/(2*n);
QoI_h = N;


end