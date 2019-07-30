clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../B-A/')

warning('off', 'stats:kmeans:EmptyCluster')

nvars = 20;
bc = 'fixed';

pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.3;
cfn = @(n) constfn(n, damping);

measParam = 0.1;

preprocfn = @(data) data;

% Create connectivity matrix.
prob = 0.5;
strength = 0.1;
force = 50;

deltat = 0.1;

numPerts = nvars;

corrThresh = 0.2;
method = 'corr';

% Perform inference
numMats = 100;
accHist = nan([numPerts, numMats]);
wccs = nan([nvars, numMats]);
predMatHists = nan(nvars, nvars, numPerts, numMats);
%hubNodes = [10];
for k = 1 : numMats
    k
    while true
        %seed = MakeNetworkSymmER(floor(nvars / 4), 0.5, true);
        %symmMat = SFNG(nvars, 1, seed);
        %mat = MakeDirected(symmMat, 1);
        mat = MakeNetworkER(nvars, 0.5, true);
        %mat(1 : 8, 1 : 8) = MakeNetworkER(8, 1, true);
        %mat(6 : 10, 6 : 10) = MakeNetworkER(5, 1, true);
        %mat(:, hubNodes) = 1;
        %mat(1 : nvars + 1 : nvars^2) = 0;

        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength * K;

        G = digraph(mat.');
        wcc = centrality(G, 'outdegree');
        wccs(:, k) = wcc;
        [sortedWCC, sortedIdx] = sort(wcc, 'descend');
        pertIdx = sortedIdx.'; %randsample(1 : nvars, nvars);

        % If this adjacency matrix is bad, make a new simulation.
        [disconnectedNodes, amplitudes, waitTime] = checkHarmonicMat(K, damping, force);
        if waitTime > 500 || ~isempty(disconnectedNodes) || any(amplitudes > -0.00001)
            continue
        end
        
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
        data = GenerateHarmonicData(nvars, tSpan, 1, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
        noisyData = WhiteGaussianNoise(data, measParam);
        
        obsIdx = true([1, nvars]);
        leftPad = 100;
        rightPad = pertLength;
        [~, predMatHist, ~, ~, ~, tableResults] = ...
            PerturbationBaseExperiment(noisyData, mat, 1, preprocfn, ...
                obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, corrThresh);
        
        predMatHists(:, :, :, k) = predMatHist;
        accHist(:, k) = sum(sum((predMatHist == repmat(mat, [1, 1, numPerts])) .* repmat(~eye(nvars), [1, 1, numPerts]), 1), 2) / (nvars^2 - nvars);
        break
    end
end

meanAccHist = mean(accHist, 2);

figure(1)
plot(1 : numPerts, meanAccHist)

figure(2)
diffAccs = diff(meanAccHist);
plot(2 : numPerts, diffAccs)
%degreeDist = prob * ones([1, numPerts - 1]);
%degreeDist(hubNodes - 1) = (nvars - 1) / nvars;
%hold on; plot(2 : numPerts, degreeDist);


plot(1 : numPerts, accRandom, 'LineWidth', 3);
hold on;
plot(1 : numPerts, accCloseness, 'LineWidth', 3);
hold on;
plot(1 : numPerts, accDegree, 'LineWidth', 3)
legend({'Random', 'Closeness', 'Degree'}, 'Fontsize', 20)
ylim([0.6 1])
