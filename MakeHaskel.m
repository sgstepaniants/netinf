function hu = MakeHaskel(u, s)
    % u(:, :, trial) = [[x_1(1), ..., x_1(m)]; ...; [x_n(1), ..., x_n(m)]]
    % hu(:, :, trial) = [[x_1(1), ..., x_1(s), ..., x_n(1), ..., x_n(s)];
    %                   ...; [x_1(m-s+1), ..., x_1(m), ..., x_n(m-s+1), ..., x_n(m)]].'

    % stack data matrix into Haskel time series matrix
    [n, m, N] = size(u);
    hu = zeros(n * s, m - s + 1, N);

    for k = 1 : m - s + 1
        u_s = u(:, k : k + s - 1, :);
        perm = permute(u_s, [2, 1, 3]);
        hu(:, k, :) = reshape(perm, [], 1, N);
    end
end
