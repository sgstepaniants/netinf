function pertValues = PertMeanVariances(data, pertTimes, movvarWidth, leftPad, rightPad, thresh)
    [nvars, timeLength] = size(data);
    numPerts = length(pertTimes);
    
    pertValues = zeros(numPerts, nvars);
    for k=1:numPerts
        leftBound = max(pertTimes(k) - leftPad, 1);
        rightBound = min(pertTimes(k) + rightPad, timeLength);
        windowData = data(:, leftBound : rightBound);
        windowDataMovvar = movvar(windowData, [movvarWidth, 0], 0, 2);
        
        meanVars = mean(windowDataMovvar, 2);
        meanVars = meanVars / max(meanVars);
        values = 1 ./ meanVars;
        values(meanVars < thresh) = Inf;
        pertValues(k, :) = values;
    end
end
