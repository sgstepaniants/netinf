function p = randfn(n, lo, hi)
    p = (hi - lo) * rand([n, 1]) + lo;
end
