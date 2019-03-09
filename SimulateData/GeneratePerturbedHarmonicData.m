function [data, pertTimes, pertLength] = GeneratePerturbedHarmonicData(mat, K, deltat, N, pfn, vfn, mfn, damping, bc, pertIdx, pertForce)
%
% calls nncoupled.m to return solution of NNCoupled model
%
% INPUTS:
%
% N
%       a scalar value for the number of random trials requested 
%
%
% K
%       a matrix of spring coefficents (the underlying connectivity matrix)
%
% pfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial positions, one for each node in the 
%       network
%
% vfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial velocities, one for each node in the 
%       network
%
% mfn 
%       a function that takes a scalar n and retuns an [n x 1] vector 
%       of (possibly random) masses, one for each node
%       in the network
%
%
% cfn 
%       a function that takes a scalar n and retuns a vector of length
%       n - 1, n, or n + 1 of (possibly random) damping coefficients,
%       one for each node in the network
%
% bc
%       a string that specifies the boundary conditions of the nearest
%       neighbor spring model ('circ', 'fixed', 'free')
%
% pertForce
%       a number that signifies the force with which the nodes are
%       perturbed in the network
%
% OUTPUTS:
%
% data
%       an [n x m x N] matrix of the solution of the ODE for n nodes,
%       m time points, and N random trials
    nvars = size(mat, 1);
    numPerts = length(pertIdx);

    % Check if this adjacency matrix has disconnected oscillators.
    G = digraph(mat.');
    distLeft = distances(G, 1);
    distRight = distances(G, nvars);
    disconnectedNodes = find(~isfinite(distLeft) & ~isfinite(distRight));

    % Check if this adjacency matrix is resonant.
    A = K(2:nvars+1, 2:nvars+1);
    A = A - diag(sum(K(2:nvars+1, :), 2));
    lambdas = eig(A);
    amplitudes = real(-damping + sqrt(damping^2 + 4 * lambdas)) / 2;

    % If this adjacency matrix is bad, make a new simulation.
    if ~isempty(disconnectedNodes)
        error('Network has disconnection components')
    end
    
    if any(amplitudes > -0.00001)
        error('Linear system has nonnegative eigenvalues')
    end

    % Create the timespan for the simulation.
    eps = 0.01;
    waitTime = ceil(log(eps / min(sqrt(sum((pertForce * inv(A)).^2)))) / max(amplitudes));
    if waitTime > 500
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
