%%% Using the modified collocation
clear all

% Define problem dimension
d = 2;

% Define number of output parameters of fvarout
OutputLength = 4;
Depth = 7; % collocation level, starting from 1


% Create full grid for plotting
gs = 33;
[X,Y] = meshgrid(linspace(0,2,gs),linspace(-1,1,gs));
XX = X(:);
YY = Y(:);
% fexact = 1./((X-(1:nout)/nout).^2+(Y).^2+(1:nout)./nout);


% Get options structure for sparse interpolation
options = spset('RelTol', 0,'AbsTol', 0,'DimensionAdaptive', 'off',...
    'Vectorized', 'off','functionArgType','vec',...
    'keepGrid','off','keepFunctionValues', 'off',...
    'MaxDepth', Depth-1,'NumberOfOutputs', 1,'OutputLength',OutputLength,...
    'VariablePositions',[1:d],'SparseIndices','off',...
    'gridtype','Maximum','enableDCT', 'off'); %clenshaw-curtis,Maximum,NoBoundary,chebyshev,Gauss-Patterson
% Here note that NumberOfOutputs = 1 always, since we treat one vector as
% an item.

% Compute sparse grid weights over domain [0,2]x[-1,1]
z = spvals(@fvarout, d, [0 2; -1 1], options,OutputLength);


% Compute inpterpolated values at full grid
ip = spinterp(z, [XX,YY]);


    % Plot interpolated results
for k = 1:OutputLength
	subplot(1,OutputLength,k);
	mesh(X, Y, reshape(ip(k,:)',size(X)));
	title(['interpolated out',num2str(k)]);
end






%------------------------------------------------------------------
function out = fvarout(x,y,nout)
% FVARARGOUT    function with multiple output arguments, the number
% of output arguments is determined by the parameter nout.
% out: vector of size nout x 1.

out = zeros(1,nout);
out(1:nout) = 1./((x-(1:nout)/nout).^2+(y).^2+(1:nout)./nout);

end