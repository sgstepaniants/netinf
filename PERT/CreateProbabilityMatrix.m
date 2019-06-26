function [AprobHist, pertOrders] = CreateProbabilityMatrix(observedData, pertIdx, obsIdx, pertTimes, ...
                                                        leftPad, rightPad, method, thresh, movvarWidth)
    % Check that all the necessary input arguments are specified.
    if nargin < 7
        error('Not enough input arguments.')
    end
    
    if nargin < 8
        movvarWidth = 0;
    end
    
    nvars = length(obsIdx);
    
    % Get the number of perturbations and the number
    % of observations.
    numPerts = length(pertIdx);
    
    % Find causal connections from perturbation data.
    % pertOrders will have dimension (numObs x nvars).
    pertOrders = GetPertOrders(observedData, pertIdx, obsIdx, pertTimes, leftPad, rightPad, method, thresh, movvarWidth);
    
    % Create a penalty factor that will multiplicatively decrease
    % probabilities of non-causal edges.
    penalty = 0.1;
    
    % Prior probability that any edge exists is 0.5
    initMat = 0.5 * ones(nvars, nvars);
    initMat(1 : nvars + 1 : nvars^2) = NaN;
    AprobHist = repmat(initMat, [1, 1, numPerts + 1]);
    AprobHist(~obsIdx, :, :) = NaN;
    AprobHist(:, ~obsIdx, :) = NaN;
    obsInds = find(obsIdx);
    for numPert=1:numPerts
        pertInd = pertIdx(numPert);
        
        % If the initially perturbed node in the cascade cannot be
        % observed, we can still find the nodes that it causes.
        if ~ismember(pertInd, find(obsIdx))
            AprobHist(pertOrders(numPert, :) == 2, pertInd, numPert + 1) = 1;
        end
        
        for j=obsInds
            for k=obsInds
                AprobHist(k, j, numPert + 1) = AprobHist(k, j, numPert);
                
                % j is the index causing and k is the index being caused.
                jPertOrder = pertOrders(numPert, j);
                kPertOrder = pertOrders(numPert, k);
                
                % If j never becomes perturbed, move on.
                if isinf(jPertOrder)
                    continue
                end
                
                causalInds = pertOrders(numPert, :) == kPertOrder - 1;
                if ~isinf(kPertOrder) && causalInds(j)
                    % If j causes k, then recompute the probability of the
                    % edge j to k.
                    causalProbs = AprobHist(k, causalInds, numPert);
                    AprobHist(k, j, numPert + 1) = AprobHist(k, j, numPert) / (1 - prod(1 - causalProbs));
                else
                    % If we have evidence that j was perturbed and k did
                    % not become perturbed immediately after, then penalize this edge.
                    if (jPertOrder < kPertOrder - 1)
                        AprobHist(k, j, numPert + 1) = penalty * AprobHist(k, j, numPert);
                    end
                end
            end
        end
    end
    
    AprobHist = AprobHist(:, :, 2 : end);
end
