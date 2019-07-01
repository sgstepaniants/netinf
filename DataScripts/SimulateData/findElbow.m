function elbowInd = findElbow(x, y)
    numpts = length(x);
    p1 = [x(1), y(1), 0];
    p2 = [x(end), y(end), 0];
    proj = cross(repmat(p2 - p1, [numpts, 1]), [x.', y.', zeros(numpts, 1)] - repmat(p1, [numpts, 1]), 2) / norm(p2 - p1);
    dists = sqrt(sum(proj.^2, 2));
    [~, elbowInd] = max(dists);
end

