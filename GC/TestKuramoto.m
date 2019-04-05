% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;
addpath('./InitFunctions/')

nvars = 3;

endtime = 10;
nobs = 200;
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);

A = MakeNetworkER(nvars, 0.5, true); % adjacency matrix
K = 1; % connection strength

pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]
cfn = @(n) constcfn(n, 1);

kickTime = 50;
pert = 1;
Y = GenerateKuramotoData(A, tSpan, 1, K, pfn, wfn, cfn, kickTime, pert)
%X = noisefn(Y);

for t = 1 : size(tSpan, 2)
    for i = 1 : nvars
        pos = 2 * pi * Y(i, t, 1);
        r = 0.1;
        rectangle('Position', [cos(pos) - r, sin(pos) - r, 2 * r, 2 * r], ...
            'FaceColor',[i / nvars, 0.2, 0.2], 'Curvature', [1, 1])
        axis([-2 2 -2 2])
    end
    pause(0.001)
    clf
end
