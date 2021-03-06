clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../B-A/')

warning('off', 'stats:kmeans:EmptyCluster')

expNum = 'PertOrders';

nvars = 10;
bc = 'fixed';

pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.3;
cfn = @(n) constfn(n, damping);

measParam = 0.1;

preprocfn = @(data) downsample(data.', 20).';

% Create connectivity matrix.
prob = 0.5;
strength = 0.1;
force = 50;

deltat = 0.1;

numPerts = nvars;

rhoThresh = 0.995;

numMats = 1;
numTypes = 3;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);
if exist(expPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
    if upper(m) == 'N'
        return
    end
    rmdir(expPath, 's')
end
mkdir(expPath)

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

tprLog = nan([1, numPerts * numTypes * numMats]);
fprLog = nan([1, numPerts * numTypes * numMats]);
accLog = nan([1, numPerts * numTypes * numMats]);
predMats = cell(1, numPerts * numTypes * numMats);
diagnosticsLog = nan([numPerts * numTypes * numMats, 3]);
%hubNodes = [10];

parDataSave = @(fname, noisyData, mat, K)...
            save(fname, 'noisyData', 'mat', 'K');    
parResultsSave = @(fname, est, tpr, fpr, acc, diagnostics)...
            save(fname, 'est', 'tpr', 'fpr', 'acc', 'diagnostics');

M = 12;
c = progress(numPerts * numTypes * numMats);
parfor (idx = 1 : numPerts * numTypes * numMats, M)
    [k, t, m] = ind2sub([numPerts, numTypes, numMats], idx);
    fprintf('perts: %d, type: %d, mat: %d\n', k, t, m)
    
    % Count the number of iterations done by the parfor loop
    c.count();
    
    currExpPath = sprintf('%s/pert%d/type%d/mat%d', expPath, k, t, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    end
    
    while true
        %seed = MakeNetworkSymmER(floor(nvars / 4), 0.5, true);
        %symmMat = SFNG(nvars, 1, seed);
        %mat = MakeDirected(symmMat, 1);
        mat = MakeNetworkER(nvars, 0.5, true);
        %mat(1 : 4, 1 : 4) = MakeNetworkER(4, 1, true);
        %mat(:, hubNodes) = 1;
        %mat(1 : nvars + 1 : nvars^2) = 0;

        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength * K;
        
        pertIdx = randsample(1 : nvars, k);
        if t > 1
            G = digraph(mat.');
            if t == 2
                wcc = centrality(G, 'outcloseness');
            elseif t == 3
                wcc = centrality(G, 'outdegree');
            end
            [~, sortedIdx] = sort(wcc, 'descend');
            pertIdx = sortedIdx.';
        end
        
        % If this adjacency matrix is bad, make a new simulation.
        [disconnectedNodes, amplitudes, waitTime] = checkHarmonicMat(K, damping, force);
        if waitTime > 500 || ~isempty(disconnectedNodes) || any(amplitudes > -0.00001)
            continue
        end
        
        endtime = waitTime * (k + 1);
        nobs = round(endtime / deltat);
        tSpan = linspace(0, endtime, nobs);
        
        % Build up forcing function.
        times = round(linspace(0, nobs, numPerts+2));
        pertTimes = times(2:end-1);
        pertLength = round(waitTime / (4 * deltat));

        forcingFunc = zeros([nvars, nobs]);
        for p=1:k
            forcingFunc(pertIdx(p), pertTimes(p):pertTimes(p)+pertLength) = force;
        end

        % Generate data with forced perturbations.
        data = GenerateHarmonicData(nvars, tSpan, 1, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
        noisyData = WhiteGaussianNoise(data, measParam);
        
        dataObsIdx = true([1, nvars]);
        [est, tableResults] = GrangerBaseExperiment(data, ...
                mat, preprocfn, dataObsIdx, rhoThresh);
        if isnan(est)
            continue
        end
        
        parDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, mat, K);
        parResultsSave(sprintf('%s/results.mat', currExpPath), est, ...
            tableResults.tpr, tableResults.fpr, tableResults.acc, tableResults.diagnostics);
        
        predMats{idx} = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

predMats = reshape(predMats, [numPerts, numTypes, numMats]);
tprLog = reshape(tprLog, [numPerts, numTypes, numMats]);
fprLog = reshape(fprLog, [numPerts, numTypes, numMats]);
accLog = reshape(accLog, [numPerts, numTypes, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numPerts, numTypes, numMats, 3]);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', ...
    'accLog', 'diagnosticsLog');


%% Plot Results

accRandom = nanmean(squeeze(accLog(:, 1, :)), 2);
accCloseness = nanmean(squeeze(accLog(:, 2, :)), 2);
accDegree = nanmean(squeeze(accLog(:, 3, :)), 2);

plot(1 : numPerts, accRandom, 'LineWidth', 3);
hold on;
plot(1 : numPerts, accCloseness, 'LineWidth', 3);
hold on;
plot(1 : numPerts, accDegree, 'LineWidth', 3)
legend({'Random', 'Closeness', 'Degree'}, 'Fontsize', 20)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
xlim([1, 20])
ylim([0.5, 1])


h = plot(G, 'LineWidth', 2, 'ArrowSize', 10);
h.Marker = 's';
p.MarkerSize = 7;
h.NodeColor = 'r';
xd = get(h, 'XData');
yd = get(h, 'YData');
nl = h.NodeLabel;
h.NodeLabel = '';
text(xd, yd, nl, 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
