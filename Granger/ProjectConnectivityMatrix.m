function G = ProjectConnectivityMatrix(U, S, G_hat)
% Take connectivity matrix of hankel matrix V time series and transform it
% into a connectivity matrix of the original X time series using weights
% from the U and S matrices where X = USV'.

[n, ~, N, reps] = size(U);

% create connectivity matrix for original X Hankel time series data
G = zeros(n, n, N, reps);
for r = 1 : reps
    for trial = 1 : N
        %G(:, :, trial, r) = U(:, :, trial, r) * G_hat(:, :, r) * U(:, :, trial, r)';
        G(:, :, trial, r) = U(:, :, trial, r) * S(:, :, trial, r) * G_hat(:, :, r) * inv(S(:, :, trial, r)) * U(:, :, trial, r)';
    end
end

end
