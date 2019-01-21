function p = unifpfn(n, bc)
    if strcmp(bc, 'circ')
        p = (0 : n-1)' / n - 0.5;
    elseif strcmp(bc, 'fixed')
        p = (1 : n)' / (n+1) - 0.5;
    elseif strcmp(bc, 'free')
        p = (0 : n-1)' / (n-1) - 0.5;
    else
        error('invalid boundary conditions')
    end
end
