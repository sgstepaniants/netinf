function pertValues = PertCorrs(data, pertIdx, pertTimes, leftPad, rightPad, corrThresh)
    [nvars, timeLength] = size(data);
    numPerts = length(pertTimes);
    
    pertValues = zeros(numPerts, nvars);
    for k=1:numPerts
        pertInd = pertIdx(k);
        
        leftBound = max(pertTimes(k) - leftPad, 1);
        rightBound = min(pertTimes(k) + rightPad, timeLength);
        pertWindow = leftBound:rightBound;
        windowData = data(:, pertWindow);
        
        corrMat = abs(corr(windowData.'));
        corrs = corrMat(pertInd, :);
        values = 1 ./ corrs;
        values(corrs < corrThresh) = Inf;
        pertValues(k, :) = values;
    end
end
