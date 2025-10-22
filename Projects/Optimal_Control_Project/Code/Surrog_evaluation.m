d=5;
ColDepth = 2;
noise_level = 0.02;
ref_para = [-0.04, 0.8, 0.012, -0.01, 0.6];
delta_para = abs(ref_para)*noise_level;
range = [ref_para - delta_para, ref_para + delta_para];
QoI_l = 602; % x & y coordinate
SaveData = false;



NoRandPts = 50;
sample = [repmat(ref_para,NoRandPts,1).*((rand(NoRandPts,5)-0.5)*2*noise_level+1)]'; %  5 x NoRandPts
% Evaluate_surrogate = spinterp(z,sample);
% figure(1),hold on; for i = 1:NoRandPts, plot(Evaluate_surrogate(1:301,i), Evaluate_surrogate(302:end,i),'r-'), end

for i=1:NoRandPts

    Evaluate_surrogate = spinterp(z,sample(:,i));
    output = direct_solve(sample(:,i)');
    figure(2),hold on;
    plot(Evaluate_surrogate(1:301), Evaluate_surrogate(302:end),'r-')
    plot(output(1:301), output(302:end),'b-.')
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


% % plot of function evaluation at sparse grids
% figure(1),hold on; for i = 1:61, plot(data(i,1:301), data(i,302:end),'r-'), end
% for i=1:3
%     NoCol = size(z.fvals{i},2);
%     for j=1:NoCol
%     plot(z.fvals{i}(1:301,j),z.fvals{i}(302:end,j),'b-.')
%     end
% end