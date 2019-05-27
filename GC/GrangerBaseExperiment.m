function [est, tableResults] = GrangerBaseExperiment(data, mats, preprocfn, dataObsIdx, rhoThresh)

    nvars = size(mats, 1); % number of variables / oscillators
    numMats = size(mats,3); % number of matrices we try
    if nargin < 4
        dataObsIdx = true([numMats, nvars]);
    end
    if nargin < 5
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

        if ~isnan(est(obsIdx, obsIdx, j))
            tableResultsDiagnostics(j, :) = diagnostics;
            tableResultsTPR(j) = nnz((est(:, :, j) + truth == 2) .* ~eye(nvars)) / numPositives;
            tableResultsFPR(j) = nnz((est(:, :, j) - truth == 1) .* ~eye(nvars)) / numNegatives;
            tableResultsAcc(j) = nnz((est(:, :, j) == truth) .* ~eye(nvars)) / (numObs^2-numObs);
        end

        count = count + 1;
    end
    
    tableResults = struct('tpr', tableResultsTPR, 'fpr', tableResultsFPR, ...
        'acc', tableResultsAcc, 'diagnostics', tableResultsDiagnostics);
end
