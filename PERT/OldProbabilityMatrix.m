function [Aprob, pertOrders] = ChngptProbabilityMatrix(data, pertIdx, obsIdx, pertTimes, changepointThresh)
    n = size(data, 1);
    numPerts = length(pertIdx);
    
    % Find causal connections from perturbation data.
    pad = 100;
    pertOrders = ChngptPertOrders(data, pertIdx, obsIdx, pertTimes, pad, changepointThresh);
    
    Aprob = ones(n) / 2;
    penalty = 0.1;
    for numPert=1:numPerts
        for j=1:n
            for k=1:n
                % k is the index causing and j is the index being caused.
                jPertTime = pertOrders(numPert, j);
                kPertTime = pertOrders(numPert, k);
                
                % If j or k cannot be observed or k never becomes perturbed, move on.
                if isnan(jPertTime) || isnan(kPertTime) || isinf(kPertTime)
                    continue
                end
                
                causalInds = pertOrders(numPert, :) == jPertTime - 1;
                if ~isinf(jPertTime) && causalInds(k)
                    % If k causes j, then recompute the probability of the
                    % edge k to j.
                    causalProbs = Aprob(j, causalInds);
                    Aprob(j, k) = Aprob(j, k) / (1 - prod(1 - causalProbs));
                else
                    % If we have evidence that k was perturbed and j did
                    % not immediately become perturbed, then penalize this edge.
                    if (kPertTime < jPertTime - 1)
                        Aprob(j, k) = penalty * Aprob(j, k);
                    end
                end
            end
        end
    end
end
