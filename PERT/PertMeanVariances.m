function pertValues = PertMeanVariances(data, pertTimes, pertLength, movvarWidth, pad, thresh)
    n = size(data, 1);
    numPerts = length(pertTimes);
    
    pertValues = zeros(numPerts, n);
    for k=1:numPerts
        leftBound = max(pertTimes(k) - pad, 1);
        rightBound = pertTimes(k) + pertLength;
        pertWindow = leftBound:rightBound;
        windowData = data(:, pertWindow);
        windowDataMovvar = movvar(windowData, [movvarWidth, 0], 0, 2);
        
        meanVars = mean(windowDataMovvar, 2);
        values = 1 ./ meanVars;
        values(meanVars < thresh) = Inf;
        pertValues(k, :) = values;
    end
end
