% MCMC for Inverse Problem in MATLAB
% Estimating parameters of a linear model with noise

clear all;
close all;
clc;

% Generate synthetic data
trueParams = [2.5, 1.2]; % True parameters [slope, intercept]
sigma = 0.5;             % Noise standard deviation

x = linspace(0, 5, 30)';
y_true = trueParams(1) * x + trueParams(2);
y_obs = y_true + sigma * randn(size(x));

% Plot true model and observed data
figure;
plot(x, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True Model');
hold on;
plot(x, y_obs, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Observed Data');
xlabel('x');
ylabel('y');
legend('Location', 'northwest');
title('True Model and Observed Data');
grid on;

% MCMC parameters
nIterations = 50000;  % Number of MCMC iterations
burnIn = 2000;        % Burn-in period
thin = 5;             % Thinning interval

% Prior distribution parameters (normal priors)
priorMean = [0, 0];   % Mean for [slope, intercept]
priorStd = [5, 5];    % Standard deviation for [slope, intercept]

% Initialize parameters
currentParams = [0, 0]; % Start with initial guess

% Store chains
chain = zeros(nIterations, 2);
acceptance = zeros(nIterations, 1);

% Calculate current likelihood
currentLikelihood = computeLikelihood(currentParams, x, y_obs, sigma);

% MCMC loop
for i = 1:nIterations
    % Propose new parameters (random walk proposal)
    proposedParams = currentParams + 0.1 * randn(1, 2);
    
    % Calculate likelihood for proposed parameters
    proposedLikelihood = computeLikelihood(proposedParams, x, y_obs, sigma);
    
    % Calculate prior probabilities
    priorCurrent = normpdf(currentParams(1), priorMean(1), priorStd(1)) * ...
                   normpdf(currentParams(2), priorMean(2), priorStd(2));
    priorProposed = normpdf(proposedParams(1), priorMean(1), priorStd(1)) * ...
                    normpdf(proposedParams(2), priorMean(2), priorStd(2));
    
    % Calculate acceptance probability
    acceptanceRatio = (proposedLikelihood * priorProposed) / (currentLikelihood * priorCurrent);
    acceptanceProb = min(1, acceptanceRatio);
    
    % Decide whether to accept proposal
    if rand() < acceptanceProb
        currentParams = proposedParams;
        currentLikelihood = proposedLikelihood;
        acceptance(i) = 1;
    end
    
    % Store current values
    chain(i, :) = currentParams;
    
    % Display progress every 1000 iterations
    if mod(i, 1000) == 0
        fprintf('Iteration: %d/%d\n', i, nIterations);
    end
end

% Calculate acceptance rate
acceptanceRate = mean(acceptance(burnIn:end));
fprintf('Acceptance rate: %.2f%%\n', acceptanceRate * 100);

% Remove burn-in and thin the chain
chain = chain(burnIn:thin:end, :);
slopeChain = chain(:, 1);
interceptChain = chain(:, 2);

% Parameter estimates
slopeEst = mean(slopeChain);
interceptEst = mean(interceptChain);
slopeStd = std(slopeChain);
interceptStd = std(interceptChain);

fprintf('Parameter estimates:\n');
fprintf('Slope = %.3f ± %.3f (true = %.3f)\n', slopeEst, slopeStd, trueParams(1));
fprintf('Intercept = %.3f ± %.3f (true = %.3f)\n', interceptEst, interceptStd, trueParams(2));

% Plot MCMC chains
figure;
subplot(2, 2, 1);
plot(slopeChain, 'b');
ylabel('Slope');
title('MCMC Chain for Slope Parameter');
grid on;

subplot(2, 2, 3);
plot(interceptChain, 'r');
ylabel('Intercept');
xlabel('Iteration');
title('MCMC Chain for Intercept Parameter');
grid on;

% Plot posterior distributions
subplot(2, 2, 2);
histogram(slopeChain, 50, 'FaceColor', 'b', 'EdgeColor', 'none', 'Normalization', 'pdf');
hold on;
xline(trueParams(1), 'r', 'LineWidth', 2, 'DisplayName', 'True Value');
xline(slopeEst, 'k--', 'LineWidth', 2, 'DisplayName', 'Estimated Value');
xlabel('Slope');
ylabel('Probability Density');
title('Posterior Distribution of Slope');
legend;

subplot(2, 2, 4);
histogram(interceptChain, 50, 'FaceColor', 'r', 'EdgeColor', 'none', 'Normalization', 'pdf');
hold on;
xline(trueParams(2), 'r', 'LineWidth', 2, 'DisplayName', 'True Value');
xline(interceptEst, 'k--', 'LineWidth', 2, 'DisplayName', 'Estimated Value');
xlabel('Intercept');
ylabel('Probability Density');
title('Posterior Distribution of Intercept');
legend;

% Plot fitted model with uncertainty
figure;
plot(x, y_obs, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Observed Data');
hold on;
plot(x, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True Model');

% Draw samples from posterior to show uncertainty
nSamples = 100;
sampleIndices = randperm(length(slopeChain), nSamples);
for i = 1:nSamples
    idx = sampleIndices(i);
    ySample = slopeChain(idx) * x + interceptChain(idx);
    s = plot(x, ySample, 'b-', 'HandleVisibility', 'off');
    s.Color(4) = 0.05;
end


% Plot estimated model
yEst = slopeEst * x + interceptEst;
plot(x, yEst, 'b--', 'LineWidth', 2, 'DisplayName', 'Estimated Model');

xlabel('x');
ylabel('y');
legend('Location', 'northwest');
title('Model Fit with Uncertainty');
grid on;

% Likelihood function
function likelihood = computeLikelihood(params, x, yObs, sigma)
    yPred = params(1) * x + params(2);
    residuals = yObs - yPred;
    likelihood = exp(-0.5 * sum((residuals / sigma).^2));
end