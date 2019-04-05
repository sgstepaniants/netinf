% Simulate data of two sine curves with (possibly) different amplitudes,
% phases, and frequencies. Fit VAR models to this data to see if they
% capture the periodicity of these functions.

clear all; close all; clc;

% Create sine data.
fn1 = @(x) sin(x);

measParam = 0.001;
T = 200;
tSpan = linspace(0, 4*pi, T);
Y = fn1(tSpan);
%X = AddGaussianNoise(Y, measParam);
X = Y;

% Compute the AIC and BIC for this model order.
%p = floor(T/2);
%[aic,bic,moaic,mobic] = tsdata_to_infocrit(X, p, 'LWR');


% Fit VAR model to data.
p = floor(T/2);
H = zeros(T-p, p);
for k=1:size(H, 1)
    H(k, :) = X((T-k):-1:(T-k-p+1));
end
y = X(T:-1:p+1).';
lasso_A = zeros(1, 1, p);
B = lasso(H, y);
lasso_A(1, 1, :) = B(:, 1);

res = y - H * squeeze(lasso_A);
covar = res.' * res / (T-p-1);


figure(2)
[pred, pred_peek] = ComputeVarFit(X, lasso_A);
plot(X.')
hold on
plot(pred.', 'o')
plot(pred_peek.', 'x')
