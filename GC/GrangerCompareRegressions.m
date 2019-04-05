% This script samples two oscillatory curves where one is a shited version
% of the other. It then fits VAR models to both these curves using OLS and
% LWR using the MVGC toolbox. It also fits VAR models using LASSO regression
% which (apart from OLS and LWR) records periodicity information from the
% data in its fit coefficients. Finally, we compare the Granger causuality
% weights in both directions for the two curves.

fn1 = @(x) sin(x);
fn2 = @(x) sin(x-pi/4);

measParam = 0.05;
T = 800;
tSpan = linspace(0, 10*pi, T);

Y = [fn1(tSpan); fn2(tSpan)];
X = AddGaussianNoise(Y, measParam);
n = size(X, 1);

figure(1)
plot(X.')

% Which curve should we predict using data from the other?
% fn2 starts oscillating later than fn1 which implies that fn1 causes fn2.
% Let's use the data from fn1 to predict fn2.
pred_idx = 2;

p = 80;

% Fit VAR models using OLS and find the optimal model order using AIC and
% BIC.
[ols_aic, ols_bic, ols_moaic, ols_mobic] = tsdata_to_infocrit(X, p, 'OLS');
tsdata_to_infocrit(X, ols_moaic, 'OLS');
load ols_A.mat
ols_A = A(pred_idx, :).';
load ols_E.mat
ols_E = E(pred_idx, :);
ols_MSE = sum(ols_E.^2)/length(ols_E)

ols_H = zeros(T-ols_moaic, n*ols_moaic);
for k=1:size(ols_H, 1)
    ols_H(k, :) = reshape(X(:, (T-k):-1:(T-k-ols_moaic+1)).', [], 1);
end
ols_y = X(pred_idx, T:-1:(ols_moaic+1)).';
ols_y_pred = ols_H * ols_A;

% Plot how OLS regression predicted the pred_idx curve.
figure(2)
plot(flip(ols_y))
hold on
plot(flip(ols_y_pred), 'x')


% Fit VAR models using LWR and find the optimal model order using AIC and
% BIC.
[lwr_aic, lwr_bic, lwr_moaic, lwr_mobic] = tsdata_to_infocrit(X, p, 'LWR');
tsdata_to_infocrit(X, lwr_moaic, 'LWR');
load lwr_A.mat
lwr_A = A(1, :).';
load lwr_E.mat
lwr_E = E(pred_idx, :);
lwr_MSE = sum(lwr_E.^2)/length(lwr_E)

lwr_H = zeros(T-lwr_moaic, n*lwr_moaic);
for k=1:size(lwr_H, 1)
    lwr_H(k, :) = reshape(X(:, (T-k):-1:(T-k-lwr_moaic+1)).', [], 1);
end
lwr_y = X(pred_idx, T:-1:(lwr_moaic+1)).';
lwr_y_pred = lwr_H * lwr_A;

% Plot how LWR regression predicted the pred_idx curve.
figure(3)
plot(flip(lwr_y))
hold on
plot(flip(lwr_y_pred), 'x')


% Fit a VAR model using LASSO regression.
H = zeros(T-p, n*p);
for k=1:size(H, 1)
    H(k, :) = reshape(X(:, (T-k):-1:(T-k-p+1)).', [], 1);
end
y = X(pred_idx, T:-1:(p+1)).';
[B, FitInfo] = lasso(H, y, 'CV', 10);
lasso_A = B(:, round(0.9*size(B, 2)));

y_pred = H * lasso_A;
lasso_E = y - y_pred;
lasso_MSE = sum(lasso_E.^2)/length(lasso_E)

% Plot how LASSO regression predicted the pred_idx curve.
figure(4)
plot(flip(y))
hold on
plot(flip(y_pred), 'x')
