function pertValues = PertCorrs(data, pertIdx, pertTimes, pertLength, pad, corrThresh)
    n = size(data, 1);
    numPerts = length(pertTimes);
    
    pertValues = zeros(numPerts, n);
    for k=1:numPerts
        pertInd = pertIdx(k);
        
        leftBound = max(pertTimes(k) - pad, 1);
        rightBound = pertTimes(k) + pertLength; %round(pertLength / 2);
        pertWindow = leftBound:rightBound;
        windowData = data(:, pertWindow);
        
        corrMat = abs(corr(windowData.'));
        corrs = corrMat(pertInd, :);
        values = 1 ./ corrs;
        values(corrs < corrThresh) = Inf;
        pertValues(k, :) = values;
    end
end
