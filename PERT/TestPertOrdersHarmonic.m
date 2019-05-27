% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../kmeans_opt/')

warning('off', 'stats:kmeans:MissingDataRemoved')

nvars = 20;
bc = 'fixed';

endtime = 200;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.1);

pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.2;
cfn = @(n) constfn(n, damping);


prob = 0.5;
spring = 0.1;
damping = 0.3;


% Number of simulation trials.
numTrials = 100;
% Number of perturbations per simulation trial.
numPerts = 1;

corrThresh = 0.5;

% Save true perturbation cascade orders.
truePertOrdersLog = zeros([numTrials, nvars]);
% Save predicted perturbation cascade orders.
predPertOrdersLog = zeros([numTrials, nvars]);

% Number of correctly predicted elements in each perturbation cascade.
correctPredElem = 0;
% Number of correctly predicted perturbation cascades.
correctPredSim = 0;

trial = 1;
while trial <= numTrials
    trial
    
    % Create connectivity matrix.
    mat = MakeNetworkER(nvars, prob, true);
    
    K = MakeNetworkTriDiag(nvars + 2, false);
    K(2:nvars+1, 2:nvars+1) = mat;
    K = spring * K;
    
    % Check if this adjacency matrix has disconnected oscillators.
    G = digraph(mat.');
    distLeft = distances(G, 1);
    distRight = distances(G, nvars);
    disconnectedNodes = find(~isfinite(distLeft) & ~isfinite(distRight));

    % Check if this adjacency matrix is resonant.
    A = K(2:nvars+1, 2:nvars+1);
    A = A - diag(sum(K(2:nvars+1, :), 2));
    lambdas = eig(A);
    amplitudes = real(-damping + sqrt(damping^2 + 4 * lambdas)) / 2;

    % If this adjacency matrix is bad, make a new simulation.
    if ~isempty(disconnectedNodes) || any(amplitudes > 0.00001)
        continue
    end

    % Create forcing function.
    forcingFunc = zeros([nvars, length(tSpan)]);
    pertIdx = randsample(nvars, numPerts);
    times = round(linspace(0, length(tSpan), numPerts+2));
    pertTimes = times(2:end-1);
    pertLength = round(length(tSpan)/(3*(numPerts+1)));
    pertForce = 30;
    for k=1:numPerts
        forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
    end

    % Simulate harmonic oscillator movement.
    ntrials = 1;
    data = GenerateHarmonicData(nvars, tSpan, numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);

    % Compute perturbation order and compare to true perturbation order.
    pad = 100;
    obsIdx = true([1, nvars]);
    observedData = data(obsIdx, :);
    predPertOrders = GetPertOrders(observedData, nvars, pertIdx, obsIdx, pertTimes, 'corr', corrThresh, pad);
    truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
    
    predPertOrdersLog(trial, :) = predPertOrders;
    truePertOrdersLog(trial, :) = truePertOrders;
    
    predPertOrders(isnan(predPertOrders)) = -1;
    truePertOrders(isnan(truePertOrders)) = -1;
    predPertOrders(isinf(predPertOrders)) = -2;
    truePertOrders(isinf(truePertOrders)) = -2;
    
    correctPredElem = correctPredElem + nvars - nnz(predPertOrders - truePertOrders);
    correctPredSim = correctPredSim + ~any(predPertOrders - truePertOrders);
    
    trial = trial + 1;
end

accuracyElem = correctPredElem / (nvars * numPerts * numTrials)
accuracySim = correctPredSim / numTrials


%% Get the TPR for every distance.

% Save TPR, FPR, FNR, TNR for every distance away from the perturbed node.
distanceStats = nan([4, nvars]);

for n=1:nvars
    if ismember(n, truePertOrdersLog)
        % Compute TPR
        distanceStats(1, n) = nnz(predPertOrdersLog(truePertOrdersLog == n) == n) / nnz(truePertOrdersLog == n);
        % Compute FPR
        distanceStats(2, n) = nnz(predPertOrdersLog(truePertOrdersLog ~= n) == n) / nnz(truePertOrdersLog ~= n);
        % Compute FNR
        distanceStats(3, n) = 1 - distanceStats(1, n);
        % Compute TNR
        distanceStats(4, n) = 1 - distanceStats(2, n);
    end
end

figure(1)
plot(distanceStats(1:2, :).')
title(sprintf('Performance of Change Point Detection for Varying Distances (n=%d)', nvars))
xlabel('Distance From Perturbed Node')
ylabel('TPR / FPR')
legend(['TPR'; 'FPR'])
