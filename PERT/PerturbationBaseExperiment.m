function [predMats, predMatsHist, AprobHist, truePertOrders, predPertOrders, tableResults] ...
    = PerturbationBaseExperiment(data, mats, numTrials, preprocfn, dataObsIdx, ...
    dataPertIdx, dataPertTimes, leftPad, rightPad, method, thresh, movvarWidth)

    if nargin < 12
        movvarWidth = 0;
    end
    
    nvars = size(mats, 1); % number of variables / oscillators
    numMats = size(mats, 3); % number of matrices we try
    numPerts = size(dataPertIdx, 2); % number of perturbations

    % Tables to hold results
    tableResultsTPR = nan(1, numMats);
    tableResultsFPR = nan(1, numMats);
    tableResultsAcc = nan(1, numMats);
    
    % est holds predicted networks
    AprobHist = nan(nvars, nvars, numPerts, numMats, numTrials);
    predMats = nan(nvars, nvars, numMats);
    truePertOrders = nan(numPerts, nvars, numMats, numTrials);
    predPertOrders = nan(numPerts, nvars, numMats, numTrials);
    
    % Loop over the networks
    for j = 1 : numMats
        obsIdx = dataObsIdx(j, :);
        numObs = nnz(obsIdx);
        pertIdx = dataPertIdx(j, :);
        pertTimes = dataPertTimes(j, :);
        
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
            
            % Get the true perturbation orders.
            truePertOrders(:, :, j, trial) = TruePertOrders(truth, pertIdx, obsIdx);

            % Run network inference on this data
            try
                [Aprob, predPertOrders(:, :, j, trial)] = CreateProbabilityMatrix(X, ...
                                                                pertIdx, obsIdx, pertTimes, ...
                                                                leftPad, rightPad, method, ...
                                                                thresh, movvarWidth);
                
                AprobHist(:, :, :, j, trial) = Aprob;
            catch e
                fprintf('%s\n', e.identifier)
                fprintf('%s\n', e.message)
                continue
            end
        end
        
        % Average connectivity probabilites to get predicted adjacency matrix
        AprobFinal = mean(squeeze(AprobHist(:, :, end, j, :)), 3);
        predMat = double(AprobFinal > 0.5);
        predMat(isnan(AprobFinal)) = NaN;
        predMats(:, :, j) = predMat;
        
        % Save results
        tableResultsTPR(j) = nnz((predMat + truth == 2) .* ~eye(nvars)) / numPositives;
        tableResultsFPR(j) = nnz((predMat - truth == 1) .* ~eye(nvars)) / numNegatives;
        tableResultsAcc(j) = nnz((predMat == truth) .* ~eye(nvars)) / (numObs^2-numObs);
    end
    
    % Get the network reconstruction our algorithm produces.
    AprobAveHist = mean(AprobHist, 5);
    predMatsHist = double(AprobAveHist > 0.5);
    predMatsHist(isnan(AprobAveHist)) = NaN;
    
    tableResults = struct('tpr', tableResultsTPR, 'fpr', tableResultsFPR, 'acc', tableResultsAcc);
end
