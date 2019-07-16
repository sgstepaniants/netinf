function [disconnectedNodes, amplitudes, waitTime] = checkHarmonicMat(K, damping, pertForce)
    nvars = size(K, 1) - 2;
    A = K(2:nvars+1, 2:nvars+1);
    mat = A ~= 0;
    
    % Check if this adjacency matrix has disconnected oscillators.
    G = digraph(mat.');
    distLeft = distances(G, 1);
    distRight = distances(G, nvars);
    disconnectedNodes = find(~isfinite(distLeft) & ~isfinite(distRight));

    % Check if this adjacency matrix has eigenvalues larger than the threshold.
    A = A - diag(sum(K(2:nvars+1, :), 2));
    lambdas = eig(A);
    amplitudes = real(-damping + sqrt(damping^2 + 4 * lambdas)) / 2;
    
    % Create the timespan for the simulation.
    waitTime = 0;
    eps = 0.01;
    if nargin == 3
        waitTime = 2 * ceil(log(eps / min(sqrt(sum((pertForce .* pinv(A)).^2)))) / max(amplitudes));
    end
end
