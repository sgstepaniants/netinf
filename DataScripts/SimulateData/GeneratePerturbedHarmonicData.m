function [data, pertTimes, pertLength] = GeneratePerturbedHarmonicData(nvars, K, deltat, N, pfn, vfn, mfn, damping, bc, pertIdx, pertForce)
    numPerts = length(pertIdx);
    
    [disconnectedNodes, amplitudes, waitTime] = checkHarmonicMat(K, damping, pertForce);
    if ~isempty(disconnectedNodes)
        error('Network has disconnection components')
    elseif any(amplitudes > -0.00001)
        error('Linear system has nonnegative eigenvalues')
    elseif waitTime > 500
        error('Simulation timespan is too long due to insufficient damping')
    end

    endtime = waitTime * (numPerts + 1);
    nobs = round(endtime / deltat);
    tSpan = linspace(0, endtime, nobs);

    % Build up forcing function.
    times = round(linspace(0, nobs, numPerts+2));
    pertTimes = times(2:end-1);
    pertLength = round(nobs/(10*(numPerts+1)));

    forcingFunc = zeros([nvars, nobs]);
    for k=1:numPerts
        forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
    end

    % Generate data with forced perturbations.
    cfn = @(n) constfn(n, damping);
    data = GenerateNNCoupledData(nvars, tSpan, N, K, pfn, vfn, ...
        mfn, cfn, bc, forcingFunc);
end
