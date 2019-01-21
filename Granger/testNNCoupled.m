% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('./InitFunctions/')

nvars = 2;
bc = 'fixed';

endtime = 10;
nobs = 500;
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);

%pfn = @(n) unifpfn(n, bc);
pfn = @(n) unifpfn(n, bc) + [1; 0];
vfn = @(n) zeros([n, 1]);


prob = 1;
spring = 1;
%mat = MakeNetworkSymmER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars + 2, false);
%K(2, 5) = 1;
%K(2:nvars+1, 2:nvars+1) = mat;
K = spring * K;

forcingFunc = zeros([nvars, length(tSpan)]);
%forcingFunc(ceil(nvars/2), 1) = 10;

ntrials = 100;
Y = GenerateNNCoupledData(nvars, tSpan, ntrials, K, pfn, vfn, ...
        @(n)constmfn(n, 1), @(n)constcfn(n, 0), bc, forcingFunc);
%X = noisefn(Y);

return
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
            axis([-1 1 -1 1])
        end
    end
    pause(0.001)
    clf
end
