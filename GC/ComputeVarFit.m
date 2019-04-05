function [ts_var, ts_var_peek, res] = ComputeVarFit(ts_data, varMats, p)
    [n, endtime] = size(ts_data);
    
    % VAR model prediction from p initial timepoints.
    ts_var = zeros(n, endtime);
    ts_var(:, 1 : p) = ts_data(:, 1 : p);
    
    % VAR model prediction if allowed to see all real data.
    ts_var_peek = zeros(n, endtime);
    ts_var_peek(:, 1 : p) = ts_data(:, 1 : p);
    
    for t = p + 1 : endtime
        for k = 1 : p
            ts_var(:, t) = ts_var(:, t) + varMats(:, :, k) * ts_var(:, t - k);
            ts_var_peek(:, t) = ts_var_peek(:, t) + varMats(:, :, k) * ts_data(:, t - k);
        end
    end
    
    ts_var(:, 1 : p) = nan;
    ts_var_peek(:, 1 : p) = nan;
    
    % Compute the residuals for the predicted data.
    res = ts_data(p+1:end) - ts_var_peek(p+1:end);
end
