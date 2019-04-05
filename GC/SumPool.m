function mat_pool = SumPool(mat, s, stride)
    % mean pool the matrix mat with a pooling filter of size s x s
    [n, m] = size(mat);
    if n < s || m < s
        error('matrix dimensions are smaller than dimensions of pooling filter')
    end

    mat_pool = zeros((n - s) / stride, (m - s) / stride);

    for j = 0 : (n - s) / stride
        for k = 0 : (m - s) / stride
            mat_pool(j + 1, k + 1) = sum(sum(mat(stride * j + 1 : stride * j + s, stride * k + 1 : stride * k + s)));
        end
    end
end
