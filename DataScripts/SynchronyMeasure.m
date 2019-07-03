function r = SynchronyMeasure(data)
    nvars = size(data, 1);
    r = abs(sum(exp(1i * data), 1)) / nvars;
end
