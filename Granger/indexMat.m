function index = indexMat(mat)
    nvars = size(mat, 1);
    
    powerMat = reshape(cumsum(reshape(ones(nvars) - eye(nvars), 1, nvars^2)) - 1, nvars, nvars).' .* ~eye(nvars);
    index = sum(sum(mat .* (2*ones(nvars)).^powerMat .* ~eye(nvars))) + 1;
end
