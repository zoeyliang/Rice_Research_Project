function [Mesh_low,QoI,h_low,z]=load_lfm(k,interpolate)


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
z = load(fullfile('./Toy_Surrogates/',filename));


end