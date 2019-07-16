% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

nvars = 3;
bc = 'fixed';

endtime = 25;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

mfn = @(n) constfn(n, 1);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

% Specify the damping constant.
damping = 0.2;
cfn = @(n) constfn(n, damping);

prob = 0.5;
strength = 1;
mat = [[0, 0, 0]; [1, 0, 0]; [0.5, 1, 0]]; %MakeNetworkSymmER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars + 2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = strength * K;

forcingFunc = zeros([nvars, length(tSpan)]);
forcingFunc(1, 100:150) = 10;

ntrials = 1;
data = GenerateHarmonicData(nvars, tSpan, ntrials, K, pfn, vfn, ...
        mfn, cfn, bc, forcingFunc);
noisyData = noisefn(data);

plot(noisyData(:, :, 1).')

preprocfn = @(data) standardize(data);
dataObsIdx = true([1, nvars]);
rhoThresh = 1;
[est, tableResults] = GrangerBaseExperiment(noisyData, mat, preprocfn, dataObsIdx, rhoThresh)

% for t = 1 : size(tSpan, 2)
%     for i = 1 : nvars
%         if (strcmp(bc, 'circ'))
%             pos = 2 * pi * Y(i, t, 1);
%             r = 0.1;
%             rectangle('Position', [cos(pos) - r, sin(pos) - r, 2 * r, 2 * r], ...
%                 'FaceColor',[i / nvars, 0.2, 0.2], 'Curvature', [1, 1])
%             axis([-2 2 -2 2])
%         else
%             pos = Y(i, t, 1);
%             r = 0.1;
%             rectangle('Position', [pos - r, - r, 2 * r, 2 * r], ...
%                 'FaceColor',[i / nvars, 0.2, 0.2], 'Curvature', [1, 1])
%             axis([-1 1 -1 1])
%         end
%     end
%     pause(0.001)
%     clf
% end
