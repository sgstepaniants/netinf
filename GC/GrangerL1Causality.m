% This script accepts time series data of X and Y stochastic processes
% and performs Granger Causality tests in either direction using L1
% regression.

fn1 = @(x) sin(x);
fn2 = @(x) sin(x-pi/4);

measParam = 0.05;
T = 800;
tSpan = linspace(0, 10*pi, T);

Y = [fn1(tSpan); fn2(tSpan)];
Y = Y - repmat(mean(Y, 2), [1, T]);
X = AddGaussianNoise(Y, measParam);
nvars = size(X, 1);

ntrials = 10;

figure(1)
plot(X.')


% Which curve should we predict using data from the other?
% fn2 starts oscillating later than fn1 which implies that fn1 causes fn2.
% Let's use the data from fn1 to predict fn2.
from_idx = 1;
to_idx = 2;


% Fit a VAR model using LASSO regression.
p = 50;

% H containts time series data from all variables in the system
H = zeros(T-p, nvars*p);
% H_prime containts time series data from all variables except the one we
% are testing for causality (at index from_idx).
H_prime = zeros(T-p, (nvars-1)*p);

X_prime = X;
X_prime(from_idx, :) = [];
for k=1:size(H, 1)
    H(k, :) = reshape(X(:, (T-k):-1:(T-k-p+1)).', [], 1);
    H_prime(k, :) = reshape(X_prime(:, (T-k):-1:(T-k-p+1)).', [], 1);
end

y = X(to_idx, T:-1:(p+1)).';
[B, FitInfo] = lasso(H, y, 'CV', 10);
A = B(:, FitInfo.Index1SE);

figure(2)
plot(A)

[B, FitInfo] = lasso(H_prime, y, 'CV', 10);
A_prime = B(:, FitInfo.Index1SE);

figure(3)
plot(A_prime)

y_pred = H * A;
E = y - y_pred;
MSE = sum(E.^2)/length(E)

y_pred_prime = H_prime * A_prime;
E_prime = y - y_pred_prime;
MSE_prime = sum(E_prime.^2)/length(E_prime)

% Plot how LASSO regression predicted the pred_idx curve.
figure(4)
plot(flip(y))
hold on
plot(flip(y_pred), 'go')
plot(flip(y_pred_prime), 'ro')

% Compute the G causality between the from_idx variable and the to_idx
% variable.
gc = log(cov(E_prime) / cov(E))


%% Perform a Statistical Test on the Granger Causality Coefficients
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
morder = p;

F = [[NaN, -0.0015]; [0.0223, NaN]];
pval = mvgc_pval(F, morder, T, ntrials, 1, 1, nvars-2, tstat);
sig  = significance(pval, alpha, mhtc);
