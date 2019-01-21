function pertOrders = MeanVarPertOrders(data, pertIdx, pertTimes, pertLength, movvarWidth, pad, thresh)
    [n, m] = size(data);
    numPerts = length(pertTimes);
    
    pertOrders = zeros(numPerts, n);
    for k=1:numPerts
        pertInd = pertIdx(k);
        
        leftBound = max(pertTimes(k) - pad, 1);
        rightBound = pertTimes(k) + pertLength;
        pertWindow = leftBound:rightBound;
        windowData = data(:, pertWindow);
        windowDataMovvar = movvar(windowData, [movvarWidth, 0], 0, 2);
        
        pertMeans = mean(windowDataMovvar, 2);
        pertAve = pertMeans(pertInd);
        pertMeans(pertInd) = [];
        pertMeans = pertMeans / pertAve;
        
        if length(pertMeans) > 1
            [idx, C] = kmeans_opt(1 ./ pertMeans);
        else
            idx = 1;
            C = pertMeans;
        end
        
        [~, sortedIdx] = sort(C);
        [~, sortedIdx] = sort(sortedIdx);
        
        order = sortedIdx(idx).' + 1;
        order(pertMeans < thresh) = Inf;
        
        pertOrder = [order(1:pertInd-1), 1, order(pertInd:end)];
        pertOrders(k, :) = pertOrder;
    end
end
