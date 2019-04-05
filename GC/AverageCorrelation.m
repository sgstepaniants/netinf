% Take average of correlation adjacency matrices for each trial of the
% nncoupled simulation and plot it. Take correlation adjacency matrix
% of average simulation over all trials and plot it.

clear all; close all; clc;

nvars = 50;

deltat = 0.1;
endtime = 10;
nobs = endtime / deltat;
tSpan = linspace(0, endtime, nobs);

ntrials = 500;

randpfn = @(n) rand(n, 1) - 0.5; % random [-0.5, 0.5]
%randpfn = @(n) (0 : n - 1)' / n - (n - 1) / (2 * n) + randn(n, 1) / (2 * n);
randvfn = @(n) zeros(n, 1); % masses start at rest (0 velocity)
randmfn = @(n) ones(n, 1); % all masses are equal
randkfn = @(n) ones(n, 1); % all spring constants are equal
randcfn = @(n) zeros(n, 1); % undamped oscillations
bc = 'circ';

% simulate nncoupled data
X = GenerateNNCoupledData(nvars, tSpan, ntrials, randpfn, randvfn, randmfn, randkfn, randcfn, bc, 0);

% compute a correlation matrix for each trial and average them
corr_sum = zeros(nvars);
for t = 1 : ntrials
    corr_sum = corr_sum + abs(corrcoef(X(:, :, t)'));
end
ave_corr = corr_sum / ntrials;

figure
subplot(1, 2, 1)
plotmat(ave_corr)
title('Average Pearson Correlation Matrix of all Time Series Trials')
axis image


% compute the average time series simulation and get its
% correlation matrix
X_ave = mean(X, 3);
corr_ave = corrcoef(X_ave');

subplot(1, 2, 2)
mat = abs(corr_ave);
plotmat(mat)
title('Pearson Correlation Matrix of Average Time Series')
axis image
