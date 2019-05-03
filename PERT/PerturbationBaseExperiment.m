function [est, tableResults] = PerturbationBaseExperiment(data, ...
    mats, numTrials, preprocfn, dataObsIdx, dataPertIdx, dataPertTimes, ...
    dataPertLength, method, corrThresh, pad, movvarWidth, freq, resultPath)
    saveData = true;
    if nargin < 14
        saveData = false;
    end
    
    % Make directory to hold results files if one does not already exist
    if saveData && exist(resultPath, 'dir') ~= 7
        error('Result path not found: %s', resultPath)
    end
    
    nvars = size(mats, 1); % number of variables / oscillators
    numMats = size(mats,3); % number of matrices we try

    % Tables to hold results
    tableResultsTPR = nan(numTrials, numMats);
    tableResultsFPR = nan(numTrials, numMats);
    tableResultsAcc = nan(numTrials, numMats);
    
    % est holds GC's estimate of the networks
    est = nan(nvars, nvars, numTrials, numMats);

    count = 1; % number of times we have run network inference method (so know how often to save work)

    % Loop over the networks
    for j = 1 : numMats
        obsIdx = dataObsIdx(j, :);
        numObs = nnz(obsIdx);
        pertIdx = dataPertIdx(j, :);
        pertTimes = dataPertTimes(j, :);
        pertLength = dataPertLength(j);
        
        truth = mats(:, :, j);
        numPositives = nnz(truth(obsIdx, obsIdx) .* ~eye(numObs));
        numNegatives = nnz(~truth(obsIdx, obsIdx) .* ~eye(numObs));

        % For each rep, different random frequencies, initial conditions, and noise
        for trial = 1 : numTrials
            if iscell(data)
                X = data{j};
                X = preprocfn(X(obsIdx, :, trial));
            else
                X = preprocfn(data(obsIdx, :, trial, j));
            end

            % Run network inference on this data
            try
                [AprobHist, predPertOrders] = CreateProbabilityMatrix(X, ...
                                                                pertIdx, obsIdx, pertTimes, ...
                                                                pertLength, method, corrThresh, pad, movvarWidth);
            catch e
                fprintf('%s\n', e.identifier)
                fprintf('%s\n', e.message)
                continue
            end
            
            % Get the true perturbation orders.
            truePertOrders = TruePertOrders(truth, pertIdx, obsIdx);

            % Get the network reconstruction our algorithm produces.
            AprobFinal = AprobHist(:, :, end);
            predMat = double(AprobFinal > 0.5);
            predMat(isnan(AprobFinal)) = NaN;
            est(:, :, trial, j) = predMat;

            % Save results
            tableResultsTPR(trial, j) = nnz((predMat + truth == 2) .* ~eye(nvars)) / numPositives;
            tableResultsFPR(trial, j) = nnz((predMat - truth == 1) .* ~eye(nvars)) / numNegatives;
            tableResultsAcc(trial, j) = nnz((predMat == truth) .* ~eye(nvars)) / (numObs^2-numObs);

            if saveData && mod(count, freq) == 0
                % save necessary files
            end
            
            count = count + 1;
        end
    end
    
    tableResults = struct('tpr', tableResultsTPR, 'fpr', tableResultsFPR, 'acc', tableResultsAcc);

    % TODO: Save whole workspace (including all those tables of results)
    if saveData
    end
end
