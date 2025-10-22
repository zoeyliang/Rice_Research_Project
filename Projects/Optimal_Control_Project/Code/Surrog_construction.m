clear, close all

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
d=5;
ColDepth = 2;
noise_level = 0.02;
ref_para = [-0.04, 0.8, 0.012, -0.01, 0.6]';
delta_para = abs(ref_para)*noise_level;
range = [ref_para - delta_para, ref_para + delta_para];
QoI_l = 602; % x & y coordinate
SaveData = false;


for depth = ColDepth

    options = spset('Vectorized', 'on','functionArgType','vec',...
        'keepGrid','off','keepFunctionValues', 'on',...
        'MaxDepth', depth,'NumberOfOutputs', 1,'QoILength',QoI_l,...
        'VariablePositions',[1:d],'SparseIndices','off',...
        'gridtype','chebyshev','enableDCT', 'off');
    z = spvals(@direct_solve, d, range, options);
    if SaveData
        strname = ['Sur_Col',num2str(depth),'.mat'];
        save([path,strname],'z','-v7.3');
    end
end


function output = direct_solve(sample)
% INPUT:
% Solver: python solver
% Argument: sample
%
% OUTPUT: a row vector containing trajectory information

results = pyrunfile("OPT_Direct_Solver.py","z",x=sample);
output = [double(results.x_opt), double(results.y_opt)]'; % use double to convert from python array to a matlab column vector.


end


% NoRandPts = 2;
% sample = repmat(ref_para,NoRandPts,1).*((rand(NoRandPts,5)-0.5)*2*noise_level+1); % NoRandPts x 5
% 
% for i=1:NoRandPts
%     % res = pyrunfile("OPT_Direct_Solver.py","z",x=sample(i,:));
% 
%     output = direct_solve(sample(i,:));
% 
%     % Evaluate_surrogate = spinterp(z,rand_para(i,:));
% end


















