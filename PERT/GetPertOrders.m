function [pertOrders, pertValues] = GetPertOrders(observedData, pertIdx, obsIdx, pertTimes, ...
                                            leftPad, rightPad, method, thresh, movvarWidth)
    addpath('../kmeans_opt/')
    
    % Give default values for all optional arguments.
    if strcmp(method, 'corr')
        % For the correlation method, every node that is perturbed must be
        % observable.
        if ~all(ismember(pertIdx, find(obsIdx)))
            error('For correlation method, perturbed nodes must be observed.')
        end
        
        if nargin < 6
            thresh = 0.5;
        elseif nargin < 7
            pad = 100;
        end
    elseif strcmp(method, 'meanvar')
        if nargin < 6
            thresh = 0;
        elseif nargin < 7
            pad = 0;
        elseif nargin < 8
            movvarWidth = 2;
        end
    elseif strcmp(method, 'chngpt')
        if nargin < 6
            thresh = 10;
        elseif nargin < 7
            pad = 100;
        end
    else
        error('Incorrect method specified.')
    end
    
    nvars = length(obsIdx);
    
    % Record values for each perturbed node.
    numPerts = length(pertTimes);
    pertValues = nan(numPerts, nvars);
    if strcmp(method, 'corr')
        cumObs = cumsum(obsIdx);
        alignedPertIdx = cumObs(pertIdx);
        pertValues(:, obsIdx) = PertCorrs(observedData, alignedPertIdx, pertTimes, leftPad, rightPad, thresh);
    elseif strcmp(method, 'meanvar')
        pertValues(:, obsIdx) = PertMeanVariances(observedData, pertTimes, movvarWidth, leftPad, rightPad, thresh);
    elseif strcmp(method, 'chngpt')
        pertValues(:, obsIdx) = PertChangepoints(observedData, pertTimes, pad, thresh);
    end
    
    % For each perturbation trial, take the values (computed by some method
    % above) for each node and sort them. Use these sorted orders to order
    % the nodes from first perturbed to last perturbed.
    pertOrders = zeros(numPerts, nvars);
    for k=1:numPerts
        % Record the index of the initially perturbed node in this trial.
        pertInd = pertIdx(k);
        
        values = pertValues(k, :);
        % Remove the value of the initially perturbed node as this value
        % generally worsens the cluster results. Also remove all NaN and Inf values.
        clusterIdx = isfinite(values);
        clusterIdx(pertInd) = 0;
        clusterValues = values(clusterIdx);
        
        % Cluster the values of all nodes in this trial.
        if isempty(clusterValues)
            idx = [];
            C = [];
        else
            if length(unique(clusterValues)) == 1 || ...
                (length(clusterValues) == 2 && abs(clusterValues(1) - clusterValues(2)) < 0.1)
                % Edge case when clustering a list of too few numbers
                idx = ones(length(clusterValues), 1);
                C = clusterValues(1);
            else
                [idx, C] = kmeans_opt(clusterValues.');
                % Edge case when clustering a list of identical numbers
                %if range(C) / mean(C) < 0.01
                %    idx = ones(length(clusterValues), 1);
                %    C = clusterValues(1);
                %end
            end
        end

        [~, sortedIdx] = sort(C);
        [~, sortedIdx] = sort(sortedIdx);

        % Set the order of the initially perturbed node in this trial to 1.
        order = zeros(1, nvars);
        order(pertInd) = 1;
        order(clusterIdx) = sortedIdx(idx).' + 1;
        order(isinf(values)) = Inf;
        pertOrders(k, :) = order;
    end
end
