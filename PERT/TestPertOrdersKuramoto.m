% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./kmeans_opt/')

% Number of network nodes
nvars = 20;

% Probability of ER network edge connections
prob = 0.5;

% Gaussian noise function
noiseVar = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Damping coefficient
damping = 0;

% Initial conditions
pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]
cfn = @(n) constcfn(n, damping);

% Connection strength
K = 50;

% Perturbation force for oscillators
pertForce = 50;


% Number of simulation trials.
numTrials = 100;
% Number of perturbations per simulation trial.
numPerts = 1;

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
    
    % Create the timespan for the simulation.
    waitTime = 5;
    deltat = 0.1;
    endtime = waitTime * (numPerts + 1);
    nobs = round(endtime / deltat);
    tSpan = linspace(0, endtime, nobs);

    % Create forcing function.
    forcingFunc = zeros([nvars, length(tSpan)]);
    pertIdx = randsample(nvars, numPerts);
    times = round(linspace(0, length(tSpan), numPerts+2));
    pertTimes = times(2:end-1);
    pertLength = round(length(tSpan)/(3*(numPerts+1)));
    for k=1:numPerts
        forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
    end

    % Simulate harmonic oscillator movement.
    ntrials = 1;
    data = GenerateKuramotoData(mat, tSpan, ntrials, K, pfn, wfn, cfn, forcingFunc);

    % Compute perturbation order and compare to true perturbation order.
    obsIdx = true([1, nvars]);
    truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
    
    %predPertOrders = ChngptPertOrders(data, pertIdx, obsIdx, pertTimes, pad, 10);
    %predPertOrders = CorrPertOrders(data, pertIdx, obsIdx, pertTimes, pertLength / 20, pad, 0.5);
    movvarWidth = 2;
    predPertOrders = MeanVarPertOrders(data, pertIdx, obsIdx, pertTimes, pertLength, movvarWidth, 0);
    
    predPertOrders(isnan(predPertOrders)) = nvars + 1;
    truePertOrders(isnan(truePertOrders)) = nvars + 1;
    predPertOrders(isinf(predPertOrders)) = nvars + 2;
    truePertOrders(isinf(truePertOrders)) = nvars + 2;
    
    predPertOrdersLog(trial, :) = predPertOrders;
    truePertOrdersLog(trial, :) = truePertOrders;
    
    correctPredElem = correctPredElem + nvars - nnz(predPertOrders - truePertOrders);
    correctPredSim = correctPredSim + ~any(predPertOrders - truePertOrders);
    
    trial = trial + 1;
end

accuracyElem = correctPredElem / (nvars * numPerts * numTrials)
accuracySim = correctPredSim / numTrials


%% Get the TPR for every distance.

% Save TPR, FPR, FNR, TNR for every distance away from the perturbed node.
distanceStats = nan([4, nvars + 2]);

for n=1:nvars+2
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
labels = strsplit(num2str(1:nvars));
labels{end + 1} = 'NaN';
labels{end + 1} = 'Inf';
bar(distanceStats(1:2, :).')
set(gca,'xticklabel',labels)
title(sprintf('Performance of Change Point Detection for Varying Distances (n=%d)', nvars))
xlabel('Distance From Perturbed Node')
ylabel('TPR / FPR')
legend(['TPR'; 'FPR'])
