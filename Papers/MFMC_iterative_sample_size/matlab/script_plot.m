% =========================================================================
% This script performs a numerical comparison of several integer-valued 
% sample size allocation schemes:
%
% Naive floor -- Peherstorfer
% Modified -- ceil with floor -- Gruber
% Iterative -- ceil with floor -- L & H
%
% relative to the continuous (real-valued) optimal allocation for the MFMC 
% method under a fixed computational budget p.
% -------------------------------------------------------------------------
%         |=============================|
%         |                             |
%         |   /.  .\           /@  @\   |
% Author: |     l      =-->      L      |
%         |      o                X     |
%         |   BEFORE            AFTER   |
%         |=============================|_{*Jiaxing Liang*}
%
% -------------------------------------------------------------------------
% |¬_¬|~_~|@_@|*_*|._.|-_-|G_G|×_×|/_\|w_w|Q_Q|>_>|=_=|+_+|#_#|z_z|>_<|6_9|
% -------------------------------------------------------------------------
% Last modified: Oct-22,2025
% =========================================================================

clear all
close all

%% Load data
% Example 1
rho = [1,   9.9977e-01   9.9925e-01  9.9728e-01   9.8390e-01]; %rho_k
rho_p1 = [rho(2:end),0];%rho_{k}
C = [7.30e+01,7.0318e-03,1.4018e-03,5.0613e-04,2.6803e-04]; %cost
% sigma1 = 1.0840e-02; %standard deviation sigma
p_max=1e4; % total cost

% % Example 2
% rho = [1, 9.999882e-01, 9.999743e-01, 9.958253e-01];
% rho_p1 = [rho(2:end),0];
% C=[44.395, 6.8409e-01, 2.9937e-01, 1.9908e-04];
% % sigma1 = 0.03;
% p_max=1e4; % total cost

delta = rho.^2-rho_p1.^2; % dalta = rho_k^2-rho_{k+1}^2

%% Initialization
M=2e5;
p_range = linspace(sum(C),p_max,M)'; % the min budget >= sum(C)
f_real_value = zeros(M,1);
f_floor = zeros(M,1);
f_iterative = zeros(M,1);
f_diff_upper_bound = zeros(M,1);

Cost_real_value = zeros(M,1);
Cost_floor = zeros(M,1);
Cost_iterative = zeros(M,1);



%% For each budget p, generate the relative difference of f and cost
for i = 1:M

    p = p_range(i);
    % =====================================================================
    % --- Feasibility check ---
    % =====================================================================
    if p < sum(C)
        error('Budget p is insufficient: must satisfy p >= sum(C).');
    end

    % =====================================================================
    % --- Real-valued sample size: formula --- Peherstorfer
    % =====================================================================
    N_star = sqrt(delta./C)*p/sum(sqrt(C.*delta));

    % =====================================================================
    % --- Integer-valued sample size: Naive floor --- Peherstorfer
    % =====================================================================
    N_floor = floor(N_star); 
    
    % =====================================================================
    % --- Modified integer-valued sample size: Ceil & floor --- Gruber
    % =====================================================================
    j=1;
    NN = N_floor;
    while N_floor(1)<1
        NN(j) = 1; % modify the corresponding entries
        N_floor = floor(sqrt(delta(j+1:end)./C(j+1:end))*(p-sum(C(1:j).*NN(1:j)))...
            /sum(sqrt(C(j+1:end).*delta(j+1:end)))); % remove the cost of the sample with size 1, and resample.
        j = j + 1;
    end
    NN(j:end) = N_floor; % the remaining sample size will be calculate using the close-form solution formula.
    

    % =====================================================================
    % --- Integer-valued sample size: Ceil & floor --- L & H
    % Iterative method with modified 
    % =====================================================================
    N_iterative = zeros(size(N_star));
    N_iterative(1) = max(floor(N_star(1)),1); % ceil if real-valued sample size falls below 1

    
    for k=2:length(N_iterative)

        N_iterative(k) = floor(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_iterative(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end))));
        
        N_iterative(k) = max(N_iterative(k),1); % modify the sample size in case of 0.

    end


    % =====================================================================
    % --- Calculate f and cost ---
    % =====================================================================
    % Real-valued
    f_real_value(i) = sum(delta./N_star);
    Cost_real_value(i) = sum(C.*N_star);
    % upper bound of relative difference of f
    f_diff_upper_bound(i) = sum(delta./(N_star-1))-f_real_value(i);

    % Integer-valued, naive floor with modified
    f_floor(i) = sum(delta./NN);
    Cost_floor(i) = sum(C.*NN);

    % Integer-valued, iterative with modified
    f_iterative(i) = sum(delta./N_iterative);
    Cost_iterative(i) = sum(C.*N_iterative);

end



%==========================================================================
% Plots of relative difference of f and cost
%==========================================================================
figure(1);hold on;
h_max = max(max((f_floor-f_real_value)./p_range,(f_iterative-f_real_value)./p_range));
h_min = abs(min(min((f_floor-f_real_value)./p_range,(f_iterative-f_real_value)./p_range)));
loglog([sum(C),sum(C)],[h_min,h_max], ':k','LineWidth',3,'MarkerSize',10);
loglog(p_range,f_diff_upper_bound./p_range,'o-','LineWidth',3,'color',[.7 0 0],'MarkerSize',5);
loglog(p_range,(f_floor-f_real_value)./p_range,'-','LineWidth',3,'color',[.2 0.4 1],'MarkerSize',30);
p2=loglog(p_range,(f_iterative-f_real_value)./p_range,'-.','LineWidth',3,'color',[.9 0.7 0],'MarkerSize',30);
p2.Color(4) = 0.6;
% loglog(p_range,f_floor-f_iterative,'-.','LineWidth',3,'color',[.6 0.7 0],'MarkerSize',30);

xlabel('Prescribed budget $p$','FontSize',25,'Interpreter','latex'); ylabel('Relative f deviation','FontSize',25,'Interpreter','latex');
lgd = legend({'$p=\sum_{k}C_k$','Upper bound: $\frac{f(N_k^*-1)}{f(N_k^*)} - 1$',...
    ['Direct floor: $\frac{f(\lfloor N_k^*\rfloor)}{f(N_k^*)}-1$'],...
    ['Iterative: $\frac{f(\lfloor M_k^*\rfloor)}{f(N_k^*)}-1$']},'Fontsize',14,'Interpreter','latex','Location','northeast'); %northeast,southwest
lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on
xlim([sum(C)/1.09 p_max*1.1])
% ylim([1e-25 2])
hold off;




figure(2);hold on;
h_max = max(max((Cost_real_value-Cost_floor)./p_range,(Cost_real_value-Cost_iterative)./p_range));
h_min = abs(min(min((Cost_real_value-Cost_floor)./p_range,(Cost_real_value-Cost_iterative)./p_range)));
loglog([sum(C),sum(C)],[1e-13,h_max], ':k','LineWidth',3,'MarkerSize',10);
loglog(p_range(1:end),sum(C)*ones(size(p_range))./p_range,'o-','LineWidth',3,'color',[.7 0 0],'MarkerSize',5);
loglog(p_range,(Cost_real_value-Cost_floor)./p_range,'-','LineWidth',3,'color',[.2 0.4 1],'MarkerSize',30);
loglog(p_range,(Cost_real_value-Cost_iterative)./p_range,'-.','LineWidth',3,'color',[.9 0.7 0],'MarkerSize',30);
% loglog(p_range,Cost_iterative-Cost_floor,'-.','LineWidth',3,'color',[.6 0.7 0],'MarkerSize',30);

xlabel('Prescribed budget $p$','FontSize',25,'Interpreter','latex'); ylabel('Relative cost deviation','FontSize',25,'Interpreter','latex');
lgd = legend({'$p=\sum_{k}C_k$','Upper bound: $\frac{\sum_{k}C_k}{p}$',...
    ['Direct floor: $1-\frac{\sum_{k}C_k\lfloor N_k^*\rfloor}{p}$'],...
    ['Iterative: $1-\frac{\sum_{k}C_k\lfloor M_k^*\rfloor}{p}$']},'Fontsize',14,'Interpreter','latex','Location','southwest');
lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on
xlim([sum(C)/1.09 p_max*1.1])
ylim([1e-20 sum(C)*2])
hold off;