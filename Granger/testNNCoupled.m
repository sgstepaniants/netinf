% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;

nvars = 4;
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);
preprocfn = @(data) NoiseThenDetrend(data, noisefn);

%randpfn = @(n) [0 : n - 1]' / n; % evenly spaced
%randpfn = @(n) [1 : n]' / (n + 1); % evenly spaced
%randvfn = @(n) [0.5, zeros(1, n - 1)]'; % one mass perturbed
%randmfn = @(n) [100, 100, 100, 100]'
%randpfn = @(n) [rand(1) / 2 - 0.25, 0.25, rand(1) / 2 + 0.25, 0.75]';
%randpfn = @(n) [-0.5, -0.25, 0, 0.25]'
randcfn = @(n) 0 * ones(n, 1);

randkfn = @(n) 1 * ones(n, 1);

bc = 'circ';
Y = GenerateNNCoupledData(nvars, tSpan, 1, randpfn, randvfn, randmfn, randkfn, randcfn, bc, 0);
%X = preprocfn(Y);

hold on
plot(Y(:, :, 1)')

for t = 1 : size(tSpan, 2)
    for i = 1 : nvars
        if (strcmp(bc, 'circ'))
            pos = 2 * pi * Y(i, t, 1);
            r = 0.1;
            rectangle('Position', [cos(pos) - r, sin(pos) - r, 2 * r, 2 * r], ...
                'FaceColor',[i / nvars, 0.2, 0.2], 'Curvature', [1, 1])
            axis([-2 2 -2 2])
        else
            pos = Y(i, t, 1);
            r = 0.1;
            rectangle('Position', [pos - r, - r, 2 * r, 2 * r], ...
                'FaceColor',[i / nvars, 0.2, 0.2], 'Curvature', [1, 1])
            axis([-0.25 1.25 -0.75 0.75])
        end
    end
    pause(0.001)
    clf
end