function p = onepfn(n, pert, bc)
    p = unifpfn(n, bc);
    idx = randi(n);
    p(idx) = p(idx) + pert * randn;
end
