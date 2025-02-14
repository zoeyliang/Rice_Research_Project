function [Mesh_low,QoI,h_low,z]=load_lfm_toy(k,interpolate)
%function [Mesh_low,QoI,h_low,z]=load_lfm_toy(k,interpolate)
%INPUT:
%k:             integer. levels of low fidelity models
%interpolate:   boolean. If the surrogate has been interpolated to the
%               common mesh or not.
%
%OUTPUT:
%Mesh_low:      vector. mesh of low fidelity model.
%QoI:           scalar. # of grid nodes of Mesh_low.
%h_low:         scalar. mesh size of Mesh_low.
%z:             structure. low fidelity surrogate.
%
% Feb-14,2025



%low-didelity mesh
level =10-k;

n = 10*2^level;
N = 2*n+1; % #of mesh grids
Mesh_low = linspace(0,2,N)';
QoI = N;
h_low = 2/(2*n);

if interpolate
    filename = strcat('Surrogate_l',num2str(level),'_interp.mat'); %'_uniform.mat'
else
    filename = strcat('Surrogate_l',num2str(level),'.mat');
end
z = load(fullfile('/Users/jiaxingliang/Google Drive/PlasmaProject_result/Toy_Surrogates/',filename));


end