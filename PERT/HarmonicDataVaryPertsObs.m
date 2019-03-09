clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./BAgraph_dir/')
addpath('./kmeans_opt/')

%% Initialize Parameters

% Number of network nodes
nvars = 5;
% Boundary conditions of oscillator system
bc = 'fixed';

% Probability of ER network edge connections
prob = 0.5;

% Gaussian noise function
noiseVar = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Spring constants in oscillator network
spring = 0.1;
% Damping in oscillator network
damping = 0.3;

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Delta t
deltat = 0.1;

% Perturbation force for oscillators
pertForce = 30;

% Threshold for correlation algorithm
corrThresh = 0.5;

% Padding for window
pad = 100;

% Require condition that all perturbed nodes must be observed
pertIsObs = true;

% Number of random matrices to generate for each choice of
% numObs and numPerts
numMats = 10;

% Number of trial simulations to make with the same random matrix
numTrials = 10;

method = 'corr';

% Make directory to hold results files if one does not already exist
expName = sprintf('EXP(nvars%d_prob%.2f_spring%.2f_damping%.2f_pertf%.2f)', nvars, prob, spring, damping, pertForce);
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


%% Simulate Data with Varying Numbers of Perturbations

fprintf('Simulate Data:\n')
for numObs = 1:nvars
    numObs
    
    % Create directory for all experiments with numObs nodes observed.
    mkdir(sprintf('%s/numobs%d', expPath, numObs))
    
    % If we can only perturb the nodes we can observe, then the maximum
    % number of nodes we can perturb is the number of observed nodes.
    maxNumPerts = nvars;
    if pertIsObs
        maxNumPerts = numObs;
    end
    
    for numPerts = 1:maxNumPerts
        numPerts
        
        % Create subdirectory for all experiments with numObs nodes
        % observed and numPerts nodes perturbed.
        mkdir(sprintf('%s/numobs%d/numperts%d', expPath, numObs, numPerts))
        
        % Create data structures to save simulation data and parameters.
        dataLog = cell(1, numMats);
        trueMats = nan(numMats, nvars, nvars);
        dataObsIdx = false(numMats, nvars);
        dataPertIdx = nan(numMats, numPerts);
        dataPertTimes = nan(numMats, numPerts);
        dataPertLength = nan(numMats);
        
        matCount = 1;
        while matCount <= numMats
            % Randomly choose which nodes we are allowed to observe.
            inds = randsample(nvars, numObs);
            obsIdx = false([1, nvars]);
            obsIdx(inds) = 1;

            % If we can only perturb the nodes we can observe, choose the
            % nodes to perturb out of the nodes we can observe.
            if pertIsObs
                if nnz(obsIdx) == 1
                    % Randsample doesn't work in the edge case where obsIdx
                    % contains one element.
                    pertIdx = find(obsIdx);
                else
                    pertIdx = randsample(find(obsIdx), numPerts);
                end
            else
                pertIdx = randsample(nvars, numPerts);
            end
            
            
            % Build up network connectivity
            mat = MakeNetworkER(nvars, prob, true);
            K = MakeNetworkTriDiag(nvars + 2, false);
            K(2:nvars+1, 2:nvars+1) = mat;
            K = spring * K;
            
            try
                [data, pertTimes, pertLength] = GeneratePerturbedHarmonicData(mat, K, deltat, numTrials, pfn, vfn, mfn, damping, bc, pertIdx, pertForce);
            catch
                continue
            end
            noisyData = noisefn(data);

            dataLog{matCount} = noisyData;
            trueMats(matCount, :, :) = mat;
            dataObsIdx(matCount, :) = obsIdx;
            dataPertIdx(matCount, :) = pertIdx;
            dataPertTimes(matCount, :) = pertTimes;
            dataPertLength(matCount) = pertLength;

            matCount = matCount + 1;

            % Save experiment simulated data, connectivity matrices, and parameters.
            save(sprintf('%s/numobs%d/numperts%d/dataLog.mat', expPath, numObs, numPerts), 'dataLog');
            save(sprintf('%s/numobs%d/numperts%d/trueMats.mat', expPath, numObs, numPerts), 'trueMats');
            save(sprintf('%s/numobs%d/numperts%d/dataObsIdx.mat', expPath, numObs, numPerts), 'dataObsIdx');
            save(sprintf('%s/numobs%d/numperts%d/dataPertIdx.mat', expPath, numObs, numPerts), 'dataPertIdx');
            save(sprintf('%s/numobs%d/numperts%d/dataPertTimes.mat', expPath, numObs, numPerts), 'dataPertTimes');
            save(sprintf('%s/numobs%d/numperts%d/dataPertLength.mat', expPath, numObs, numPerts), 'dataPertLength');
        end
    end
end
