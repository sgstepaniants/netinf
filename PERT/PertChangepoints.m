function [pertOrders, pertResponseTimes] = PertChangepoints(data, pertIdx, obsIdx, pertTimes, pad, changepointThresh)
    [n, m] = size(data);
    numPerts = length(pertTimes);
    
    pertResponseTimes = zeros(numPerts, n);
    newPertTimes = [1, pertTimes, m];
    for k=2:numPerts+1
        leftBound = max(newPertTimes(k-1), newPertTimes(k)-pad);
        rightBound = newPertTimes(k+1);
        pertWindow = leftBound:rightBound;
        
        pertResponseTimes(k-1, :) = GetResponseTimes(data(:, pertWindow), newPertTimes(k)-pertWindow(1)+1, obsIdx, changepointThresh);
    end
    
    pertOrders = zeros(numPerts, n);
    for k=1:numPerts
        pertResponseTime = pertResponseTimes(k, :).';
        order = zeros(n, 1);
        order(~isfinite(pertResponseTime)) = pertResponseTime(~isfinite(pertResponseTime));
        
        % If there is only one unique non-NaN and non-Inf element in
        % pertResponseTime, then don't perform clustering.
        finiteElems = pertResponseTime(isfinite(pertResponseTime));
        uniqueElems = unique(finiteElems);
        if ~isempty(uniqueElems)
            if length(uniqueElems) == 1
                uniqueElem = uniqueElems(1);
                order(pertResponseTime == uniqueElem) = 1;
            else
                [idx, C] = kmeans_opt(finiteElems);
                [~, sortedIdx] = sort(C);
                [~, sortedIdx] = sort(sortedIdx);
                order(isfinite(pertResponseTime)) = sortedIdx(idx);
            end
        end
        
        % Make sure that the perturbed node is the only node with label 1.
        if sum(order == 1) > 1
            order(pertIdx(k)) = 0;
            order = order + 1;
        end
        
        pertOrders(k, :) = order;
    end
end
