function pertOrders = TruePertOrders(mat, pertIdx, obsIdx)
    n = size(mat, 1);
    numPerts = length(pertIdx);
    
    % Get the distances of all nodes in the graph from the perturbed nodes.
    G = digraph(mat.');
    dist = distances(G, pertIdx);
    
    % Label all nodes that we cannot observe with distance NaN from the
    % perturbed node.
    pertOrders = nan(numPerts, n);
    for k = 1:numPerts
        pertInd = pertIdx(k);
        row = dist(k, obsIdx);
        unq = unique(row);
        [~, pertOrders(k, obsIdx)] = ismember(row, unq);
        
        % If the initially perturbed node is not in the list of observed
        % nodes, label it with a 1.
        if pertOrders(k, pertInd) ~= 1
            pertOrders(k, obsIdx) = pertOrders(k, obsIdx) + 1;
            pertOrders(k, pertInd) = 1;
        end
    end
    
    % All nodes not in the same connected component as the perturbed nodes
    % in pertIdx are labeled with distance Inf.
    pertOrders(isinf(dist)) = Inf;
end
