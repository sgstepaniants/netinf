% Draw Time Series Trajectories for Paper

clear all; close all; clc;

addpath('SimulateData/')
addpath('SimulateData/InitFunctions/')
addpath('../GC/')

% Specify noise for data.
measParam = 1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

tSpan = 0:0.01:10;

f1 = @(t) sin(2*t) + cos(t).^2 - 3 * sin(t + pi/8).^4 + 7 * sin(t).^4 .* cos(t - pi/5).^2;
f2 = @(t) -2 * sin(t) + 3 * cos(t).^2 - 6 * pi * sin(t + pi/8).^4 + 7 * sin(t).^4 .* cos(t - pi/5).^2;
f3 = @(t) 4 * sin(t) + cos(t).^2 - 3 * sin(t + pi/8).^4 + 7 * sin(5*t).^4 .* cos(t - pi/5).^2;
f4 = @(t) sin(t) + 4 * cos(3*t).^2 - 3 * sin(t + pi/8).^4 + 7 * sin(t).^4 .* cos(t - pi/5).^2;
f5 = @(t) sin(t) + cos(t).^2 - 3 * sin(2*t + pi/8).^4 + 7 * pi * sin(t).^4 .* cos(t - pi/5).^2;

% Clean unperturbed data
figure(1)
plot(tSpan, normalize(f1(tSpan)) + 4, ...
    tSpan, normalize(f2(tSpan)) + 2, ...
    tSpan, normalize(f3(tSpan)), ...
    tSpan, normalize(f4(tSpan)) - 2, ...
    tSpan, normalize(f5(tSpan)) - 4, 'LineWidth', 10)
axis([0, tSpan(end), -5, 6])
set(gca, 'xtick', [], 'ytick', [])

% Noisy unperturbed data
figure(2)
plot(tSpan, noisefn(normalize(f1(tSpan))) + 4, ...
    tSpan, noisefn(normalize(f2(tSpan))) + 2, ...
    tSpan, noisefn(normalize(f3(tSpan))), ...
    tSpan, noisefn(normalize(f4(tSpan))) - 2, ...
    tSpan, noisefn(normalize(f5(tSpan))) - 4, 'LineWidth', 3)
axis([0, tSpan(end), -5, 6])
set(gca, 'xtick', [], 'ytick', [])

% Noisy perturbed data
figure(3)
plot(tSpan, noisefn(normalize(f1(tSpan) + 10 * exp(-50 * (tSpan - .5).^2) + 10 * exp(-50 * (tSpan - 5).^2) + 10 * exp(-50 * (tSpan - 8).^2))) + 4, ...
    tSpan, noisefn(normalize(f2(tSpan) + 30 * exp(-50 * (tSpan - 1).^2) + 30 * exp(-50 * (tSpan - 5.5).^2) + 30 * exp(-50 * (tSpan - 8.5).^2))) + 2, ...
    tSpan, noisefn(normalize(f3(tSpan) + 30 * exp(-50 * (tSpan - 2).^2) + 30 * exp(-50 * (tSpan - 4.5).^2) + 30 * exp(-50 * (tSpan - 7.5).^2))), ...
    tSpan, noisefn(normalize(f4(tSpan) + 20 * exp(-50 * (tSpan - 1.5).^2) + 20 * exp(-50 * (tSpan - 5).^2) + 20 * exp(-50 * (tSpan - 8).^2))) - 2, ...
    tSpan, noisefn(normalize(f5(tSpan) + 35 * exp(-50 * (tSpan - 1.5).^2) + 35 * exp(-50 * (tSpan - 4).^2) + 35 * exp(-50 * (tSpan - 9).^2))) - 4, 'LineWidth', 3)
axis([0, tSpan(end), -5, 6])
set(gca, 'xtick', [], 'ytick', [])


nvars = 5;
tSpan = 0:0.1:25;
forcingFunc = zeros(nvars, length(tSpan));
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Noisy harmonic oscillator data
bc = 'fixed';
mfn = @(n) constfn(n, 1);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

damping = 0;
cfn = @(n) constfn(n, damping);
K = [[0, 1, 0, 0, 0, 0, 0]; [1, 0, 1, 0, 0, 0, 0]; [1, 1, 0, 0, 1, 0, 0]; [0, 1, 0, 0, 1, 1, 0]; [0, 0, 0, 1, 0, 0, 1]; [0, 0, 0, 0, 0, 0, 1]; [0, 0, 0, 0, 0, 1, 0]];
K = (K + K.') > 0

data = GenerateHarmonicData(nvars, tSpan, 1, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
figure(4)
plot(noisefn(normalize(data)).' + repmat([2; 1; 0; -1; -2], [1, length(tSpan)]).', 'LineWidth', 3);

% Noisy Kuramoto oscillator data
A = [[0, 0, 1, 0, 0]; [1, 0, 0, 1, 0]; [0, 0, 0, 0, 1]; [0, 1, 1, 0, 0]; [0, 1, 0, 0, 0]];
K = 1;
wfn = @(n) randfn(n, -1, 1);

data = GenerateKuramotoData(A, tSpan, 1, K, pfn, wfn, forcingFunc);
figure(4)
plot(noisefn(normalize(data)).' + repmat([2; 1; 0; -1; -2], [1, length(tSpan)]).', 'LineWidth', 3);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca,'TickLength',[0 0])

% Plot confusion matrices
GCHar = [[100, 0, 0, 0]; [2, 98, 0, 0]; [0, 0, 99, 1]; [0, 0, 0, 100]];
figure(5)
imagesc(GCHar)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap pink

GCKur = [[71, 10, 12, 7]; [2, 3, 0, 95]; [1, 0, 2, 97]; [1, 2, 0, 97]];
figure(6)
imagesc(GCKur)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap pink

CCMHar = [[13, 7, 1, 79]; [23, 49, 1, 27]; [18, 2, 48, 32]; [7, 1, 0, 92]];
figure(7)
imagesc(CCMHar)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap pink

CCMKur = [[90, 2, 8, 0]; [7, 3, 0, 90]; [3, 2, 2, 93]; [1, 1, 2, 96]];
figure(8)
imagesc(CCMKur)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap pink
