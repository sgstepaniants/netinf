function p = uniffn(n, lo, hi, bc)
    if strcmp(bc, 'circ')
        p = (hi - lo) * (0 : n-1)' / n + lo;
    elseif strcmp(bc, 'fixed')
        p = (hi - lo) * (1 : n)' / (n+1) + lo;
    elseif strcmp(bc, 'free')
        p = (hi - lo) * (0 : n-1)' / (n-1) + lo;
    else
        error('invalid boundary conditions')
    end
end
