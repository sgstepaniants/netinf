% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../BAgraph_dir/')

nvars = 15;

endtime = 10;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn = @(data) WhiteGaussianNoise(data, 0.1);

prob = 0.5;
mat = MakeNetworkER(nvars, prob, true);
%mat(1:floor(nvars/2), 1:floor(nvars/2)) = BAgraph_dir(floor(nvars/2), 2, floor(nvars/4)).';
%mat(floor(nvars/2)+1:end, floor(nvars/2)+1:end) = BAgraph_dir(floor(nvars/2), 2, floor(nvars/4)).';
mat(2:4, 2:4) = 1;
%mat(12:14, 12:14) = 1;
mat(1:nvars+1:nvars^2) = 0;

G = digraph(mat.');
figure(1)
plot(G)

pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

strength = 100; % connection strength

numTrials = 100;
forcingFunc = zeros([nvars, nobs]);
data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);

repData = permute(repmat(data, [1, 1, 1, nvars]), [1, 4, 2, 3]);
aveCorr = mean(cos(repData - permute(repData, [2, 1, 3, 4])), 4);

figure(2)
aves = aveCorr(:, :, 100);
aves(1 : nvars + 1 : nvars^2) = [];
%h = histogram(aves, 20);
imagesc(normalize(aveCorr(:, :, 100)) > 0.85)

figure(3)
imagesc(aveCorr(:, :, 100))


%% Perturb Network in Optimal Locations

% Width of moving variance window
movvarWidth = 10;

% Threshold for mean variance algorithm
meanThresh = 0; %0.001;

method = 'meanvar';

% Preprocessing function for data.
preprocfn = @(data) data;

% Amount of time to wait between perturbations
waitTime = 10;

measParam = 0.01;

force = 100;

pertIdx = 1 : nvars;
numPerts = length(pertIdx);

deltat = 0.01;
endtime = waitTime * (numPerts + 1);
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

% Build up forcing function.
times = round(linspace(0, nobs, numPerts+2));
pertTimes = times(2:end-1);
pertLength = round(waitTime / (4 * deltat));

forcingFunc = zeros([nvars, nobs]);
for p=1:numPerts
    forcingFunc(pertIdx(p), pertTimes(p):pertTimes(p)+pertLength) = force;
end

% Generate data with forced perturbations.
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

numTrials = 100;
data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
noisyData = nan([nvars, nobs, 1, numTrials]);
noisyData(:, :, 1, :) = WhiteGaussianNoise(data, measParam);

obsIdx = true([1, nvars]);
leftPad = 0;
rightPad = pertLength;
truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
[predMats, predMatsHist, ~, ~, predPertOrders, tableResults] = ...
    PerturbationBaseExperiment(noisyData, repmat(mat, [1, 1, numTrials]), 1, preprocfn, ...
        repmat(obsIdx, [numTrials, 1]), repmat(pertIdx, [numTrials, 1]), repmat(pertTimes, [numTrials, 1]), leftPad, rightPad, method, meanThresh, movvarWidth);

accHist = squeeze(mean(sum(sum((predMatsHist == repmat(mat, [1, 1, numPerts, numTrials])) .* repmat(~eye(nvars), [1, 1, numPerts, numTrials]), 1), 2) / (nvars^2 - nvars), 4));

plot(1 : numPerts, accHist)
