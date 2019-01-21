function mat_rescaled = rescale(mat)
    mat_min = min(min(mat));
    mat_max = max(max(mat));
    mat_rescaled = (mat - mat_min) / (mat_max - mat_min);
end
