%script to build surroagte
%ColloDepth: integer.  sparse grid level
%method: character. method to build surrogate. here we use 'spgridcb'.
%d: integer. dimension of parameter space.
%perturb_ratio1: decimal. noise in parameter1.
%perturb_ratio2: decimal. noise in parameter2.
%
%

clear all;
close all


ColloDepth = 1;
d = 2;
method = 'spgridcb';
perturb_ratio1 = 0.03;
perturb_ratio2 = 0.03;
range = [[-1,1]*perturb_ratio1;[-1,1]*perturb_ratio2];

for ii = 1:8

    n = 10*2^ii;
    N = 2*n+1; % #of mesh grids
    x = linspace(0,2,N);
    h = 2/(2*n);
    QoI_l = N;

    clear u
    for k = 0:ColloDepth 
        levelseq = spgetseq(k,d);
        X{k+1} = feval(method,levelseq); % Col nodes in the parameter space w/o scaling, stored in grid cell x{k+1}
        nColPts = size(X{k+1},1); % # of collocation pts corresponding to the current level
        v=[];
        for l = 1:d
            v(:,l) = X{k+1}(:,l) .* (range(l,2) - range(l,1)) + range(l,1); % scaling col nodes
        end
        for i = 1:nColPts
            xi= v(i,1);
            eta = v(i,2);
            u(i,:) = FEM_solver(n,xi,eta);
            %         plot(x,u(i,:),'r');
            %         hold on
        end
        z{k+1} = u';
        if k>0
            ip = zeros(QoI_l,nColPts);
            ipmethod = 'spcmpvalscb';
            % find \Delta A_{q,d} (f) = \sum a*( f(x) - U^{i-1}(f)(x) )
            for m = 0:k-1
                oldlevelseq = spgetseq(m,d);
                ip = ip + feval(ipmethod, d, QoI_l, z{m+1}, X{k+1}, ...
                    levelseq, oldlevelseq);
            end
            z{k+1} = z{k+1} - ip;
        end
    end
    save(['/Users/jiaxingliang/Google Drive/PlasmaProject_result/Toy_Surrogates/Surrogate_l',num2str(ii),'.mat'],'z');
end


