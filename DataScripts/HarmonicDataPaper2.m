% Simulate unperturbed data and vary connection probabilities and coupling strengths.

clear all; close all; clc;
addpath('SimulateData/')
addpath('SimulateData/InitFunctions/')

expNum = 'Paper5';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 10;

% Initialize masses, positions, and velocities of oscillators.
mfn = @(n) randfn(n, 0.7, 1.3);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

% Specify the damping constant.
damping = 0.25;
cfn = @(n) constfn(n, damping);

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Specify forcing function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Probabilities of network connections to try.
probs = 0.05:0.05:1;
numProbs = length(probs);

% Connection strengths to try.
strengths = 0.5:0.5:10;
numStrengths = length(strengths);

% Number of matrices to try for each probability and connection strength
% combination.
numMats = 10;

% Number of simulation trials.
numTrials = 1;

% Make directory to hold data files if one does not already exist
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);
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

% Create random connectivity matrices and simulate oscillator trajectories.
dataLog = nan(nvars, nobs, numTrials, numMats, numProbs, numStrengths);
trueMats = nan(nvars, nvars, numMats, numProbs, numStrengths);
Ks = nan(nvars+2, nvars+2, numMats, numProbs, numStrengths);
for j = 1 : numProbs
    prob = probs(j)
    
    for k = 1 : numStrengths
        strength = strengths(k)
        
        l = 1;
        count = 1;
        while l <= numMats && count <= Inf
            % Create adjacency matrices.
            mat = MakeNetworkER(nvars, prob, true);
            K = MakeNetworkTriDiag(nvars+2, false);
            K(2:nvars+1, 2:nvars+1) = mat;
            K = strength .* K; %(strength + (rand(nvars+2, nvars+2)-0.5)/2) .* K;
            
            % If any nodes in the network are not connected to the walls or
            % the eigenvalues of the system have positive real parts, don't
            % use this network.
            [~, amplitudes] = checkHarmonicMat(K, damping);
            if any(amplitudes > 0) %|| ~isempty(disconnectedNodes)
                count = count + 1;
                continue
            end
            
            trueMats(:, :, l, j, k) = mat;
            Ks(:, :, l, j, k) = K;
            
            % Simulate oscillator trajectories.
            data = GenerateHarmonicData(nvars, tSpan, ...
            numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
            noisyData = noisefn(data);
            dataLog(:, :, :, l, j, k) = noisyData;
            
            l = l + 1;
        end
    end
end

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
save(sprintf('%s/Ks.mat', expPath), 'Ks');