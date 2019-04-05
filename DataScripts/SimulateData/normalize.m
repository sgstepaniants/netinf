function vv = normalize(v)
    finiteIdx = isfinite(v);
    u = v(finiteIdx);
    uu = (u - min(u)) / (max(u) - min(u));
    vv = v;
    vv(finiteIdx) = uu;
end
