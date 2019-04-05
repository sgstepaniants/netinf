clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./BAgraph_dir/')
addpath('./kmeans_opt/')

%% Initialize Parameters

% Number of network nodes
nvars = 10;
% Boundary conditions of oscillator system
bc = 'fixed';

% Probability of ER network edge connections
probs = 0.1:0.1:1;

% Gaussian noise function
noiseVar = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Spring constants in oscillator network
springs = 0.1:0.1:1;
% Damping in oscillator network
damping = 0.3;

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);
cfn = @(n) constfn(n, damping);

% Number of nodes we can observe in the network
numObs = nvars;

% Delta t
deltat = 0.1;

% Perturbation force for oscillators
pertForce = 30;

% Threshold for correlation algorithm
corrThresh = 0.5;

% Padding for window
pad = 100;

% Number of experimental trials
numTrials = 10;

method = 'corr';

% Make directory to hold results files if one does not already exist
expName = sprintf('EXPVaryStrengthsProbs(nvars%d_numObs%d_damping%.2f_pertf%.2f)', nvars, numObs, damping, pertForce);
expPath = sprintf('HarmonicExperiments/%s', expName);
if exist(expPath, 'dir') ~= 7
    mkdir(expPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
    if upper(m) == 'N'
       return
    end
end

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));


%% Siumlate Data with Varying Network Connection Probabilities and Strengths
numProbs = length(probs);
numStrengths = length(springs);

dataLog = cell(numProbs, numStrengths, numTrials);
trueMats = cell(numProbs, numStrengths, numTrials);
dataParamsLog = cell(numProbs, numStrengths, numTrials, 4);

fprintf('Simulate Data:\n')
for i = 1:numProbs
    prob = probs(i)
    
    for j = 1:numStrengths
        spring = springs(j)

        trial = 1;
        while trial <= numTrials
            % Perturb all nodes sequentially.
            pertIdx = 1:nvars;
            numPerts = length(pertIdx);
            
            % Observe numObs nodes in the network.
            inds = randsample(nvars, numObs);
            obsIdx = false([1, nvars]);
            obsIdx(inds) = 1;
            
            % Build up network connectivity
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
            if ~isempty(disconnectedNodes) || any(amplitudes > -0.00001)
                continue
            end
            
            % Create the timespan for the simulation.
            eps = 0.01;
            waitTime = ceil(log(eps / min(sqrt(sum((pertForce * inv(A)).^2)))) / max(amplitudes));
            if waitTime > 500
                continue
            end

            endtime = waitTime * (numPerts + 1);
            nobs = round(endtime / deltat);
            tSpan = linspace(0, endtime, nobs);

            % Build up forcing function.
            times = round(linspace(0, nobs, numPerts+2));
            pertTimes = times(2:end-1);
            pertLength = round(nobs/(10*(numPerts+1)));

            forcingFunc = zeros([nvars, nobs]);
            for k=1:numPerts
                forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
            end

            % Generate data with forced perturbations.
            data = GenerateNNCoupledData(nvars, tSpan, 1, K, pfn, vfn, ...
                mfn, cfn, bc, forcingFunc);
            noisyData = noisefn(data);

            dataLog{i, j, trial} = noisyData;
            trueMats{i, j, trial} = mat;
            dataParamsLog{i, j, trial, 1} = obsIdx;
            dataParamsLog{i, j, trial, 2} = pertIdx;
            dataParamsLog{i, j, trial, 3} = pertTimes;
            dataParamsLog{i, j, trial, 4} = pertLength;

            trial = trial + 1;

            % Save experiment simulated data, connectivity matrices, and parameters.
            save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
            save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
            save(sprintf('%s/dataParamsLog.mat', expPath), 'dataParamsLog');
        end
    end
end


%% Evaluate Algorithm on Data and Show Performance on Varying Numbers of Perturbations 
truePertOrdersLog = cell(numProbs, numStrengths, numTrials);
predPertOrdersLog = cell(numProbs, numStrengths, numTrials);
predMatsHist = cell(numProbs, numStrengths, numTrials);
accuracyLog = zeros(nvars, numProbs, numStrengths, numTrials);

fprintf('Run Algorithm:\n')
for j = 1:numProbs
    prob = probs(j)
    
    for k = 1:numStrengths
        spring = springs(k)
        
        trial = 1;
        while trial <= numTrials
            data = dataLog{j, k, trial};
            if isempty(data)
                trial = trial + 1;
                continue
            end
            
            mat = trueMats{j, k, trial};
            obsIdx = dataParamsLog{j, k, trial, 1};
            pertIdx = dataParamsLog{j, k, trial, 2};
            pertTimes = dataParamsLog{j, k, trial, 3};
            pertLength = dataParamsLog{j, k, trial, 4};

            % Select only the subset of nodes that we can observe.
            observedData = data(obsIdx, :);
            [AprobHist, predPertOrders] = CreateProbabilityMatrix(observedData, ...
                                                            pertIdx, obsIdx, pertTimes, ...
                                                            pertLength, method, corrThresh, pad);
            predPertOrdersLog{j, k, trial} = predPertOrders;

            % Get the true perturbation orders.
            truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
            truePertOrdersLog{j, k, trial} = truePertOrders;

            % Get the network reconstruction our algorithm produces.
            predMatHist = AprobHist > 0.5;
            predMatsHist{j, k, trial} = predMatHist;

            % Compute the accuracy of this reconstruction (fraction of edges
            % correctly predicted) and store it in the scoring matrix.
            numObs = nnz(obsIdx);
            diff = predMatHist(obsIdx, :, :) - repmat(mat(obsIdx, :), [1, 1, nvars]);
            accHist = sum(sum(diff == 0, 1), 2) / (nvars * numObs);
            accuracyLog(:, j, k, trial) = accHist;

            % Save experiment results.
            save(sprintf('%s/truePertOrdersLog.mat', expPath), 'truePertOrdersLog');
            save(sprintf('%s/predPertOrdersLog.mat', expPath), 'predPertOrdersLog');
            save(sprintf('%s/predMatsHist.mat', expPath), 'predMatsHist');
            save(sprintf('%s/accuracyLog.mat', expPath), 'accuracyLog');
            
            trial = trial + 1;
        end
    end
end

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = mean(squeeze(accuracyLog(10, :, :, :)), 3);
%clims = [0, 1];
imagesc(aveAccuracies, clims)
set(gca,'xtick',[])
set(gca,'ytick',[])
%colorbar
%set(gca, 'XTickLabel', springs)
%set(gca, 'YTickLabel', probs)
%xlabel('Spring Constants')
%ylabel('Connection Probabilities')
