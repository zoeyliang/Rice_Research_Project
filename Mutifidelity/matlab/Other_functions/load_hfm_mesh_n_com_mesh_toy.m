function [Mesh_com,Mesh_high,h_high,QoI_h] = load_hfm_mesh_n_com_mesh_toy(com_mesh_level,hfm_mesh_level)
%function [Mesh_com,Mesh_high,h_high,QoI_h] = load_hfm_mesh_n_com_mesh_toy(com_mesh_level,hfm_mesh_level)
%INPUT:
%com_mesh_level:    integer. level of common mesh.
%hfm_mesh_level:    integer. level of mesh for high fidelity model.
%
%OUTPUT:
%Mesh_com:          structure. common mesh.
%Mesh_high:         structure. mesh of high fidelity model.
%h_high:            scalar. mesh size of high fidelity model.
%QoI_h:             integer. # of mesh grids of high fidelity model.
%
% Feb-14,2025.


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