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

M=2e5;
p_range = linspace(sum(C),p_max,M)';
f_real_value = zeros(M,1);
f_floor = zeros(M,1);
f_iterative = zeros(M,1);
f_diff_upper_bound = zeros(M,1);

Cost_real_value = zeros(M,1);
Cost_floor = zeros(M,1);
Cost_iterative = zeros(M,1);

Efficiency_real_value = zeros(M,length(C));
Efficiency_floor = zeros(M,length(C));
Efficiency_iterative = zeros(M,length(C));

for i=1:M
    p=p_range(i);
    % --- Feasibility check ---
    if p < sum(C)
        error('Budget p is insufficient: must satisfy p >= sum(C).');
    end

    delta = rho.^2-rho_p1.^2; % dalta = rho_k^2-rho_{k+1}^2

    % real-valued sample size, using close-form solution
    N_star = sqrt(delta./C)*p/sum(sqrt(C.*delta));

    % =====================================================================
    % Naive floor
    N_floor = floor(N_star); 
    
    % Modified MFMC using naive floor
    j=1;
    NN = N_floor;
    while N_floor(1)<1
        NN(j) = 1; % modify the corresponding entries
        N_floor = floor(sqrt(delta(j+1:end)./C(j+1:end))*(p-sum(C(1:j).*NN(1:j)))...
            /sum(sqrt(C(j+1:end).*delta(j+1:end)))); % remove the cost of the sample with size 1, and resample.
        j=j+1;
    end
    NN(j:end) = N_floor; % the remaining sample size will be calculate using the close-form solution formula.
    % =====================================================================

    % =====================================================================
    % Iterative strategy
    N_iterative = zeros(size(N_star));
    N_iterative(1) = max(NN(1),1);

    % --- Iterative computation for k=2,...,K ---
    for k=2:length(N_iterative)

        N_iterative(k) = floor(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_iterative(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end)))); % remove the cost of the sample with size 1, and resample.

        N_iterative(k) = max(N_iterative(k),1); % modify the corresponding entries

    end
    % =====================================================================

    f_real_value(i) = sum(delta./N_star);
    Cost_real_value(i) = sum(C.*N_star);
    Efficiency_real_value(i,:) = delta./(C.*N_star.^2);

    f_floor(i) = sum(delta./NN);
    Cost_floor(i) = sum(C.*NN);
    Efficiency_floor(i,:) = delta./(C.*NN.^2);

    f_iterative(i) = sum(delta./N_iterative);
    Cost_iterative(i) = sum(C.*N_iterative);
    Efficiency_iterative(i,:) = delta./(C.*N_iterative.^2);

    f_diff_upper_bound(i) = sum(delta./(N_star-1))-f_real_value(i);

end


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
% loglog(p_range,Cost_real_value)
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



figure(3);hold on;
loglog(p_range(1:end),Efficiency_real_value(:,1),'ro')
loglog(p_range(1:end),Efficiency_floor(:,1),'bo')
loglog(p_range(1:end),Efficiency_iterative(:,1),'g-')

lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on

figure(4);hold on;
loglog(p_range(1:end),Efficiency_real_value(:,2),'ro')
loglog(p_range(1:end),Efficiency_floor(:,2),'bo')
loglog(p_range(1:end),Efficiency_iterative(:,2),'g-')
lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on

figure(5);hold on;
loglog(p_range(1:end),Efficiency_real_value(:,3),'ro')
loglog(p_range(1:end),Efficiency_floor(:,3),'bo')
loglog(p_range(1:end),Efficiency_iterative(:,3),'g-')
lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on

figure(6);hold on;
loglog(p_range(1:end),Efficiency_real_value(:,4),'ro')
loglog(p_range(1:end),Efficiency_floor(:,4),'bo')
loglog(p_range(1:end),Efficiency_iterative(:,4),'g-')
lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on

figure(7);hold on;
loglog(p_range(1:end),Efficiency_real_value(:,5),'ro')
loglog(p_range(1:end),Efficiency_floor(:,5),'bo')
loglog(p_range(1:end),Efficiency_iterative(:,5),'g-')
lgd.FontSize = 17;
set(gca,'YScale', 'log', 'XScale', 'log','FontSize',20)
grid on
box on