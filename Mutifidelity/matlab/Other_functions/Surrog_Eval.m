function ip = Surrog_Eval(z,range,d,QoI_l,rand_vec)
%function ip = Surrog_Eval(z,range,d,QoI_l,rand_vec)
%INPUT:
%z:         structure. surrogate.
%range:     dx2 vector. range of parameter space.
%d:         dimension of parameter space.
%QoI_l:     integer.length of quantity of interest.
%rand_vec:  1xd vector. samples to be evaluated.
%
%OUTPUT:
%ip:        structure.surrogate.
%
% Feb-14,2025.



ipmethod = 'spinterpcb';
nfrom = 0;
ColDepth = 1;
nn = ColDepth;



varargin = {rand_vec(1),rand_vec(2)};
% u_from_solver = FEM_solver(n,varargin{1},varargin{2});

for k = 1:length(varargin)
    if ~isempty(range)
        varargin{k} = (varargin{k}(:)-range(k,1))./(range(k,2)-range(k,1));
    else
        varargin{k} = varargin{k}(:);
    end
end
y = [varargin{:}];
ninterp = size(y,1); % number of new currents to evaluate
ip = zeros(QoI_l,ninterp);
for k = nfrom:nn
    levelseq = spgetseq(k,d);
    ip = ip + feval(ipmethod,d,QoI_l,z{k+1},y,levelseq,[]);
end






% plot(x,ip,'r');
% hold on
% plot(x,u_from_solver,'k');
%
% L2err = sqrt(sum(h*(u_from_solver - ip').^2));
% s = s + L2err;

% title('h=1.9531e-4, ColDepth = 2');
% legend('evaluation surrogate','direct solver','Location','northwest')
% L2err_m = s/M;
% L2err_m
% hold on
% plot(x,uexact2);
% find convergence rate:
% n = 10*[2,4,8,16,32,64,128];
% h = 2./(2*n);
% p = polyfit(log(h),log(relative_err),1); % p(1) is the convergence rate



%% ========================================================================