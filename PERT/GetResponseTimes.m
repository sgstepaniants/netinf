function responseTimes = GetResponseTimes(pertWindowData, pertTime, obsIdx, changepointThresh)
    n = size(pertWindowData, 1);
    
    % Compute windowed variance of data.
    windowLength = 10;
    dataMovvar = movvar(pertWindowData, [windowLength, 0], 0, 2);

    % Record response times of oscillators after perturbation.
    responseTimes = inf([1, n]);
    for k=1:n
        %cpt = pertTime - 1 + findchangepts(dataMovvar(k, pertTime:end, :), 'Statistic','linear','MinThreshold',1);
        
        cpt = findchangepts(dataMovvar(k, :, :), 'Statistic','linear','MinThreshold',changepointThresh); % harmonic
        %cpt = pertTime - 1 + find(abs(diff(dataMovvar(k, pertTime:end, :), 2)) > changepointThresh); % kuramoto
        
        if ~isempty(cpt)
            responseTimes(k) = cpt(1);
        end
    end
    
    % Label all nodes that we cannot observe with NaN.
    vec = true(1, n);
    vec(obsIdx) = 0;
    responseTimes(:, vec) = NaN;
    
    %plot(dataMovvar.'); hold on; plot(responseTimes, 0, 'r*')
end
