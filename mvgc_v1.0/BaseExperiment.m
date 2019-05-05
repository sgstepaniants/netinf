function BaseExperiment(expnum, mats, randpfn, randvfn, ...
    randmfn, randkfn, randcfn, preprocfn, deltat, endtime, ntrials, reps, ...
    tsplits, freq, stacks, networkInferenceFn, numDiagnostics)

%
% Core function of this codebase: runs a basic experiment. There are
% many inputs so that everything can be varied. Default parameters are
% saved in UsualParams.mat. See ExperimentA1.m for an example of calling
% this function.
%
% INPUTS:
%
% expnum
%       a string giving the experiment "number," such as "A1" - this 
%       determines the subdirectory where the results are saved
%
% mats
%       an [n x n x numMats] sequence of true network structures 
%       (binary adjacency matrices) to try. The outer loop in 
%       this function loops over each network to generate data
%       on that network
%
% randpfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial positions, one for each node in the 
%       network
%
% randvfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial velocities, one for each node in the 
%       network
%
% randmfn 
%       a function that takes a scalar n and retuns an [n x 1] vector 
%       of (possibly random) masses, one for each node
%       in the network
%
% randkfn 
%       a function that takes a scalar n and retuns a vector of length
%       n - 1, n, or n + 1 of (possibly random) spring constants,
%       one for each node in the network
%
% randcfn 
%       a function that takes a scalar n and retuns a vector of length
%       n - 1, n, or n + 1 of (possibly random) damping coefficients,
%       one for each node in the network
%
% preprocfn
%       a function that takes the [n x m x N] matrix of the solution Y of 
%       the ODE for n nodes, m time points, and N random trials, preprocesses
%       it, and returns X, the version of Y used as the input to the network
%       inference function. Default is that X is the noisy version of Y.
%
% deltat
%       scalar for size of time step we want in the data we generate
%
% endtime
%       scalar for the last time point to solve the system at. (solve the 
%       system of ODEs for t = [0, endtime])
%
% ntrials
%       a scalar value for the number of random trials requested 
%
% reps
%       a scalar value for the number of times to generate data and
%       apply the network inference method for a single network + K pair. 
%       We can optionally generate the data repeatedly and infer the 
%       network on each instantiation. Then the final inferred network
%       is a "vote" over these repetitions for whether or not each 
%       edge is included. Default is 1: no repetitions.
%
% tsplits
%       a vector of time indices. We can optionally split the data into
%       smaller time intervals. We infer the network on each time interval,
%       then "vote" over the results for whether or not each edge is included.
%       This vector should contain the endpoints of each interval. Default is 
%       no splitting: tpslits = nobs (number of observations / time points in
%       data).
%
% freq
%       integer for how often to save the current workspace. Default is 10: 
%       after every 10 runs of the network inference method, save the 
%       current workspace so we can check preliminary results while the
%       code continues to run. Saving the workspace also makes it easier
%       to restart the experiment close to where we left off if something
%       goes wrong with the computer.
%
% stacks
%       number of times to stack data when creating Hankel matrix
%
% networkInferenceFn
%       function that accepts the time series data as input and outputs an 
%       adjacency matrix and any number of diagnostics. 
%
% 
% numDiagnostics 
%       scalar that states how many diagnostics you expect networkInferenceFn 
%       to output.
%

if nargin <= 15
    networkInferenceFn = @(data) DemoMVGC(data);
    numDiagnostics = 3;
end

% make directory to hold results files
mkdir(sprintf('exp%s',expnum))

nvars = size(mats, 1); % number of variables / oscillators
numMats = size(mats,3); % number of matrices we try
numSplits = length(tsplits); % number of splits in time
possedges = nvars^2-nvars; % possible edges for this number of variables
numVoters = reps * numSplits; % number of estimated matrices to vote over

% tables to hold results
tableTrueEdges = zeros(numMats, 1);
tableResultsNorm = zeros(numMats, reps, numSplits); 
tableResultsNormVoting = zeros(numMats, 1);
tableResultsDiagnostics = zeros(numMats, reps, numSplits, numDiagnostics);
tableResultsInfEdges = zeros(numMats, reps, numSplits);
tableResultsInfEdgesVoting = zeros(numMats, 1);
tableResultsFalsePos = zeros(numMats, reps, numSplits);
tableResultsFalsePosVoting = zeros(numMats, 1);
tableResultsFalseNeg = zeros(numMats, reps, numSplits);
tableResultsFalseNegVoting = zeros(numMats, 1);
tableResultsPerWrong = zeros(numMats, reps, numSplits);
tableResultsPerWrongVoting = zeros(numMats, 1);

% set up time sampling
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

% hold results of all trials and repetitions
Yhist = zeros(nvars, nobs, ntrials, reps);
[stackNobs, rank, ~] = size(preprocfn(Yhist(:, :, 1, 1)));
Xhist = zeros(rank, stackNobs, ntrials, reps);
Uhist = zeros(nvars, rank, ntrials, reps);
Shist = zeros(rank, rank, ntrials, reps);
Vhist = zeros(stackNobs, rank, ntrials, reps);

% est holds network inference method's estimate of the networks
est = zeros(nvars, nvars, reps, numSplits);

count = 1; % number of times have run network inference method (so know how often to save work)

% loop over the networks
for j = 1 : numMats
    truth = mats(:, :, j);
    tableTrueEdges(j) = nnz(truth);
    
    % for each rep, different random frequencies, initial conditions, and noise
    for r = 1 : reps
        % generate the data and save information about it
        Y = GenerateNNCoupledData(nvars, tSpan, ntrials, randpfn, randvfn, randmfn, randkfn, randcfn);
        %X = preprocfn(Y);
        [U, S, V] = preprocfn(Y);
        X = permute(V, [2, 1, 3]);
        size(X)
        size(Xhist)
        % save this repetition
        Yhist(:, :, :, r) = Y;
        Xhist(:, :, :, r) = X;
        Uhist(:, :, :, r) = U;
        Shist(:, :, :, r) = S;
        Vhist(:, :, :, r) = V;
        
        % potentially split the time interval 
        beg = 1; % beginning of first time interval
        for s = 1 : numSplits
            % run network inference on this time interval
            [est(:, :, r, s), diagnostics] = networkInferenceFn(X(:, beg : tsplits(s), :));
            beg = tsplits(s) + 1; % beginning for next time interval
            count = count + 1;
            
            % save results
            tableResultsNorm(j, r, s) = norm(est(:, :, r, s) - truth) / norm(truth);
            tableResultsDiagnostics(j, r, s, :) = diagnostics;
            tableResultsInfEdges(j, r, s) = nnz(est(:, :, r, s)); 
            tableResultsFalsePos(j, r, s) = length(find(est(:, :, r, s) - truth == 1)); 
            tableResultsFalseNeg(j, r, s) = length(find(truth - est(:, :, r, s) == 1));
            tableResultsPerWrong(j, r, s) = length(find(truth ~= est(:, :, r, s))) / possedges;
            
            if mod(count,freq) == 0
                % after every freq runs of networkInferenceFn, save the current state,
                % in case want to check on preliminary results
                save(sprintf('./exp%s/exp%s-partial.mat',expnum,expnum));
            end
        end
        
        votingMat = sum(sum(est, 3), 4) / numVoters;
        tableResultsNormVoting(j) = norm(votingMat - truth) / norm(truth);
        tableResultsInfEdgesVoting(j) = nnz(votingMat);      
        tableResultsFalsePosVoting(j) = length(find(votingMat - truth == 1)); 
        tableResultsFalseNegVoting(j) = length(find(truth - votingMat == 1));
        tableResultsPerWrongVoting(j) = length(find(truth ~= votingMat)) / possedges;
        
        % save little file with results from this network
        save(sprintf('./exp%s/Exp%s-Mat%d.mat',expnum,expnum,j),'truth','X1','Y1','est','tsplits')
    end
end

% update partial results and save whole workspace (including all those tables of results)
save(sprintf('./exp%s/exp%s-partial.mat',expnum,expnum));
save(sprintf('./exp%s/exp%s.mat',expnum,expnum));
