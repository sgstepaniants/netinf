function [est, tableResults] = GrangerBaseExperiment(data, mats, preprocfn, freq, resultPath, dataObsIdx, rhoThresh)
    saveData = true;
    if nargin < 5 || strcmp(resultPath, '')
        saveData = false;
    end
    
    % The numTrials variable controls how many trials from the data you
    % actually want to analyze. With too many trials, the analysis may
    % become prohibitively slow.
    
    % Make directory to hold results files if one does not already exist
    if saveData && exist(resultPath, 'dir') ~= 7
        error('Result path not found: %s', resultPath)
    end
    
    nvars = size(mats, 1); % number of variables / oscillators
    numMats = size(mats,3); % number of matrices we try
    if nargin < 6
        dataObsIdx = true([numMats, nvars]);
    end
    if nargin < 7
        rhoThresh = 1;
    end

    % Tables to hold results
    tableResultsDiagnostics = nan(numMats, 3);
    tableResultsTPR = nan(1, numMats);
    tableResultsFPR = nan(1, numMats);
    tableResultsAcc = nan(1, numMats);
    
    % est holds GC's estimate of the networks
    est = nan(nvars, nvars, numMats);

    count = 1; % number of times we have run network inference method (so know how often to save work)

    % Loop over the networks
    for j = 1 : numMats
        obsIdx = dataObsIdx(j, :);
        numObs = nnz(obsIdx);
        
        truth = mats(:, :, j);
        numPositives = nnz(truth(obsIdx, obsIdx) .* ~eye(numObs));
        numNegatives = nnz(~truth(obsIdx, obsIdx) .* ~eye(numObs));
            
        if iscell(data)
            X = data{j};
            X = preprocfn(X(obsIdx, :, :));
        else
            X = preprocfn(data(obsIdx, :, :, j));
        end

        % Run network inference on this data
        try
            [est(obsIdx, obsIdx, j), diagnostics] = DemoMVGC(X, rhoThresh);
        catch e
            fprintf('%s\n', e.identifier)
            fprintf('%s\n', e.message)
            continue
        end

        % Save results
        if ~isnan(est(obsIdx, obsIdx, j))
            tableResultsDiagnostics(j, :) = diagnostics;
            tableResultsTPR(j) = nnz((est(:, :, j) + truth == 2) .* ~eye(nvars)) / numPositives;
            tableResultsFPR(j) = nnz((est(:, :, j) - truth == 1) .* ~eye(nvars)) / numNegatives;
            tableResultsAcc(j) = nnz((est(:, :, j) == truth) .* ~eye(nvars)) / (numObs^2-numObs);
        end

        if saveData && mod(count, freq) == 0
            % save necessary files
        end

        count = count + 1;
    end
    
    tableResults = struct('tpr', tableResultsTPR, 'fpr', tableResultsFPR, ...
        'acc', tableResultsAcc, 'diagnostics', tableResultsDiagnostics);

    % TODO: Save whole workspace (including all those tables of results)
    if saveData
    end
end
