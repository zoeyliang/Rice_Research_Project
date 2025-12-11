%% Load data
example_number = 2;


[rho, cost, budget] = MFMC_examples(example_number)

% Compute Delta_k =  (rho_{1,k}^2 - rho_{1,k+1}^2)
delta  = rho.^2 - [rho(2:end), 0].^2;  

m0 = 0;

%[x1,fval] = MFMC_real_samples_analytical(delta, cost, budget)

[x, fval] = MFMC_real_samples(delta, cost, budget, m0)
[x, fval] = MFMC_real_samples_analytical(delta, cost, budget)



%%
function [x, fval] = MFMC_real_samples(delta, cost, budget, m0)
%
% Solve real valued sample size problem
% min  sum delta(j) / m(j)
% s.t. m(1) >= m0
%      m(j) >= m(j-1), j = 1, ..., k
%      sum cost(j) *  m(j) <= budget
%
% returns solution x and corresponding objectve function value 

  k = length(delta); % number of models

  A = zeros(k+1,k);
  b = zeros(k+1,1);
  A(1,1) = -1; b(1) = -m0;              % -m(1) <= -m0
  for j = 2:k
      A(j,j-1) = 1; A(j,j) = -1;        % m(j-1) - m(j) <= 0
  end
  A(k+1,:) = cost(:)'; b(k+1) = budget; % sum m(j)*cost(j) <= budget

  options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
                         'SpecifyObjectiveGradient',true);

  x0 = ones(k,1);

  [x, fval, flag] = fmincon(@(x)obj_fun(x,delta), x0, A, b, [], [], [], [], [], options);

end

function [val, grad] = obj_fun(x, delta)
  % objective function  sum delta(j) / x(j)   and gradient
  val  = sum(delta(:)./x);
  grad = -delta(:)./(x.^2);
end


%%
function [x, fval] = MFMC_real_samples_analytical(delta, cost, budget)
%
% Analytical solution of real sample size problem
% min  sum delta(j) / m(j)
% s.t. m(1) >= 0
%      m(j) >= m(j-1), j = 1, ..., k
%      sum cost(j) *  m(j) <= budget
%
% returns solution x and corresponding objectve function value 

  k = length(delta); % number of models

  % check that conditions on correlations and costs are satisfied
  if( any( delta <= 0 ) )
      fprintf( 'delta(j) <= 0 for some j. delta = ')
      fprintf( '%10.4e,', delta); fprintf( '\n')
      return
  end

  deltac = delta./cost;
  if( any( deltac(2:k) <= deltac(1:k-1) ) )
      fprintf( 'delta(j)/cost(j) <= delta(j-1)/cost(j-1) for some j. delta/cost = ')
      fprintf( '%10.4e,', deltac);  fprintf( '\n')
      return
  end

  scal = budget/ sum(sqrt(delta.*cost));
  x    = scal * sqrt(delta./cost);
  x    = x(:);

  fval = ( sum(sqrt(delta.*cost)) )^2 / budget;

end