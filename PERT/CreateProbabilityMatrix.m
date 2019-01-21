function [Aprob, pertOrders] = CreateProbabilityMatrix(observedData, nvars, pertIdx, obsIdx, pertTimes, ...
                                                        pertLength, method, thresh, pad, movvarWidth)
    % Check that all the necessary input arguments are specified.
    if nargin < 9
        error('Not enough input arguments.')
    end
    
    if nargin < 10
        movvarWidth = 0;
    end
    
    % Get the number of perturbations and the number
    % of observations.
    numPerts = length(pertIdx);
    numObs = length(obsIdx);
    
    % Find causal connections from perturbation data.
    % pertOrders will have dimension (numObs x nvars).
    pertOrders = GetPertOrders(observedData, nvars, pertIdx, obsIdx, pertTimes, pertLength, method, thresh, pad, movvarWidth);
    
    % Initialize the probabilit matrix we will construct.
    Aprob = ones(nvars) / 2;
    % Create a penalty factor that will multiplicatively decrease
    % probabilities of non-causal edges.
    penalty = 0.1;
    
    for numPert=1:numPerts
        pertInd = pertIdx(numPert);
        
        % If the initially perturbed node in the cascade cannot be
        % observed, we can still find the nodes that it causes.
        if ~ismember(pertInd, obsIdx)
            Aprob(pertOrders(numPert, :) == 1, pertInd) = 1;
        end
        
        for j=1:numObs
            for k=1:numObs
                % j is the index causing and k is the index being caused.
                jPertOrder = pertOrders(numPert, j);
                jObsIdx = obsIdx(j);
                kPertOrder = pertOrders(numPert, k);
                kObsIdx = obsIdx(k);
                
                % If j never becomes perturbed, move on.
                if isinf(jPertOrder)
                    continue
                end
                
                causalInds = pertOrders(numPert, :) == kPertOrder - 1;
                if ~isinf(kPertOrder) && causalInds(j)
                    % If j causes k, then recompute the probability of the
                    % edge j to k.
                    causalProbs = Aprob(kObsIdx, causalInds);
                    Aprob(kObsIdx, jObsIdx) = Aprob(kObsIdx, jObsIdx) / (1 - prod(1 - causalProbs));
                else
                    % If we have evidence that j was perturbed and k did
                    % not immediately become perturbed, then penalize this edge.
                    if (jPertOrder < kPertOrder - 1)
                        Aprob(kObsIdx, jObsIdx) = penalty * Aprob(kObsIdx, jObsIdx);
                    end
                end
            end
        end
    end
end
