function [votingMats, est, tableResults] = GrangerBaseExperiment(data, mats, numTrials, reps, preprocfn, freq, resultPath, dataObsIdx)
    saveData = true;
    if nargin < 7
        saveData = false;
    end
    
    % Make directory to hold results files if one does not already exist
    if saveData && exist(resultPath, 'dir') ~= 7
        error('Result path not found: %s', resultPath)
    end
    
    nvars = size(mats, 1); % number of variables / oscillators
    numMats = size(mats,3); % number of matrices we try
    if nargin < 7
        dataObsIdx = true([numMats, nvars]);
    end

    % Tables to hold results
    tableResultsNorm = nan(reps, numMats);
    tableResultsNormVoting = nan(1, numMats);
    tableResultsDiagnostics = nan(reps, numMats, 3);
    tableResultsTPR = nan(reps, numMats);
    tableResultsTPRVoting = nan(1, numMats);
    tableResultsFPR = nan(reps, numMats);
    tableResultsFPRVoting = nan(1, numMats);
    tableResultsAcc = nan(reps, numMats);
    tableResultsAccVoting = nan(1, numMats);

    votingMats = nan(nvars, nvars, numMats);
    
    % est holds GC's estimate of the networks
    est = nan(nvars, nvars, reps, numMats);

    count = 1; % number of times we have run network inference method (so know how often to save work)

    % Loop over the networks
    for j = 1 : numMats
        obsIdx = dataObsIdx(j, :);
        numObs = nnz(obsIdx);
        
        truth = mats(:, :, j);
        numPositives = nnz(truth(obsIdx, obsIdx) .* ~eye(numObs));
        numNegatives = nnz(~truth(obsIdx, obsIdx) .* ~eye(numObs));

        % For each rep, different random frequencies, initial conditions, and noise
        for r = 1 : reps
            % Read data and split up simulation trials into repetitions
            repSize = floor(numTrials/reps);
            trialRange = (r-1)*repSize+1:r*repSize;
            
            if iscell(data)
                X = data{j};
                X = preprocfn(X(obsIdx, :, trialRange));
            else
                X = preprocfn(data(obsIdx, :, trialRange, j));
            end

            % Run network inference on this data
            try
                [est(obsIdx, obsIdx, r, j), diagnostics] = DemoMVGC(X);
            catch e
                fprintf('%s\n', e.identifier)
                fprintf('%s\n', e.message)
                continue
            end

            % Save results
            tableResultsNorm(r, j) = norm(est(:, :, r, j) - truth) / norm(truth);
            tableResultsDiagnostics(r, j, :) = diagnostics;
            tableResultsTPR(r, j) = nnz((est(:, :, r, j) + truth == 2) .* ~eye(nvars)) / numPositives;
            tableResultsFPR(r, j) = nnz((est(:, :, r, j) - truth == 1) .* ~eye(nvars)) / numNegatives;
            tableResultsAcc(r, j) = nnz((est(:, :, r, j) == truth) .* ~eye(nvars)) / (numObs^2-numObs);

            if saveData && mod(count, freq) == 0
                % save necessary files
            end
            
            count = count + 1;
        end

        votingMat = round(nanmean(est(:, :, :, j), 3));
        votingMats(:, :, j) = votingMat;
        tableResultsNormVoting(j) = norm(votingMat - truth) / norm(truth);
        tableResultsTPRVoting(j) = nnz((votingMat + truth == 2) .* ~eye(nvars)) / numPositives;
        tableResultsFPRVoting(j) = nnz((votingMat - truth == 1) .* ~eye(nvars)) / numNegatives;
        tableResultsAccVoting(j) = nnz((votingMat == truth) .* ~eye(nvars)) / (numObs^2-numObs);
    end
    
    tableResults = struct('norm', tableResultsNorm, 'normVoting', tableResultsNormVoting, ...
                    'tpr', tableResultsTPR, 'tprVoting', tableResultsTPRVoting, ...
                    'fpr', tableResultsFPR, 'fprVoting', tableResultsFPRVoting, ...
                    'acc', tableResultsAcc, 'accVoting', tableResultsAccVoting, ...
                    'diagnostics', tableResultsDiagnostics);

    % TODO: Save whole workspace (including all those tables of results)
    if saveData
    end
end
