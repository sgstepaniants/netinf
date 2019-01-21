function prob = CalculateProbability(nvars, causal_idx, causal_hist, exclude)
    prob = 0;
    
    % Remove empty cell array contents.
    causal_hist = causal_hist(~cellfun('isempty', causal_hist));
    
    num_records = length(causal_hist);
    if num_records == 0
        return
    end
    
    restrictionMat = zeros(nvars, num_records);
    for j=1:num_records
        restrictionMat(cell2mat(causal_hist(j)), j) = 1;
    end
    
    % List all possible ways to connect causal nodes to the pert_idx node.
    load combinationMat
    workingCombinationsIdx = all(combinationMat * restrictionMat, 2) & ~(combinationMat * exclude);
    
    % Calculate the probability that an edge exists.
    prob = sum(combinationMat(workingCombinationsIdx, causal_idx)) / sum(workingCombinationsIdx);
end
