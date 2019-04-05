function vv = standardize(v)
    tLen = size(v, 2);
    vv = (v - repmat(mean(v, 2), [1, tLen, 1])) ./ repmat(std(v, 0, 2), [1, tLen, 1]);
end
