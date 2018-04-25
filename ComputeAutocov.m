function autocov_mat = ComputeAutocov(mat, j)
    % calculate the jth autocovariances of the data in each row
    % (which represents a time series)
    
    if j < 0
    
    [n, m] = size(mat);
    autocov_mat = zeros(n, m);
    
    mean_mat = MovingMean(mat);
    
    sum = zeros(n);
    for k = 1 : m
        sum = sum + mat(:, k);
        mean_mat(:, k) = sum / k;
    end
end
