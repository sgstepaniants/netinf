function p = smallrandpfn(n, bc)
    p = unifpfn(n, bc);
    if strcmp(bc, 'circ')
        p = p + randn(n, 1) / (2*n);
    elseif strcmp(bc, 'fixed')
        p = p + randn(n, 1) / (2*(n+1));
    elseif strcmp(bc, 'free')
        p = p + randn(n, 1) / (2*(n-1));
    else
        error('invalid boundary conditions')
    end
end
