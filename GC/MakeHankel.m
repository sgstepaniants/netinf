function [U_trunc, S_trunc, V_trunc] = MakeHankel(X, stacks, rank)
    % X(:, :, trial) = [[x_1(1), ..., x_1(m)]; ...; [x_n(1), ..., x_n(m)]]
    % H(:, :, trial) = [[x_1(1), ..., x_1(stacks), ..., x_n(1), ..., x_n(stacks)];
    %                   ...; [x_1(m-stacks+1), ..., x_1(m), ..., x_n(m-stacks+1), ..., x_n(m)]].'
    
    % stack data matrix into Hankel time series matrix
    [n, m, N] = size(X);
    
    if nargin == 2
        rank = n;
    end
    
    H = zeros(n * stacks, m - stacks + 1, N);

    for k = 1 : m - stacks + 1
        slice = X(:, k : k + stacks - 1, :);
        perm = permute(slice, [2, 1, 3]);
        H(:, k, :) = reshape(perm, [], 1, N);
    end
    
    
    % the Hankel matrix is generally too large to process so take the SVD
    % and return the V eigenvalue time series instead
    
    U_trunc = zeros(n * stacks, rank, N);
    S_trunc = zeros(rank, rank, N);
    V_trunc = zeros(m - stacks + 1, rank, N);
    for j = 1 : N
        j
        [U, S, V] = svd(H(:, :, j), 'econ');
        U_trunc(:, :, j) = U(:, 1 : rank);
        S_trunc(:, :, j) = S(1 : rank, 1 : rank);
        V_trunc(:, :, j) = V(:, 1 : rank);
    end
end
