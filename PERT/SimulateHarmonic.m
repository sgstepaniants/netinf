% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./best_kmeans/')
addpath('./kmeans_opt/')

warning('off', 'stats:kmeans:EmptyCluster')

nvars = 20;
bc = 'fixed';

endtime = 200;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.1);

mfn = @(n) constmfn(n, 1);
pfn = @(n) randpfn(n);
vfn = @(n) zeros([n, 1]);

% Create connectivity matrix.
prob = 0.5;
spring = 1;
mat = MakeNetworkSymmER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars + 2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = spring * K;

damping = 0.3;


% Print if this adjacency matrix will give bad simulations.
G = digraph(mat.');
distLeft = distances(G, 1);
distRight = distances(G, nvars);
disconnectedNodes = find(~isfinite(distLeft) & ~isfinite(distRight))

A = K(2:nvars+1, 2:nvars+1);
A = A - diag(sum(K(2:nvars+1, :), 2));
lambdas = eig(A);
amplitudes = real(-damping + sqrt(damping^2 + 4 * lambdas)) / 2

determinant = det(A)

% Create forcing function.
forcingFunc = zeros([nvars, length(tSpan)]);
pertIdx = 2;
numPerts = length(pertIdx);
times = round(linspace(0, length(tSpan), numPerts+2));
pertTimes = times(2:end-1);
pertLength = 200;
pertForce = 30;
for k=1:numPerts
    forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
end

eps = 0.1;
minWaitTime = log(eps / min(sqrt(sum(((pertForce * inv(A)).^2))))) / max(amplitudes)
maxWaitTime = log(eps / max(sqrt(sum(((pertForce * inv(A)).^2))))) / max(amplitudes)

% Simulate harmonic oscillator movement.
ntrials = 1;
Y = GenerateNNCoupledData(nvars, tSpan, ntrials, K, pfn, vfn, ...
        mfn, @(n)constcfn(n, damping), bc, forcingFunc);

% Plot oscillator trajectories.
figure(1)
plot(Y.'); hold on; plot((-inv(A) * forcingFunc).')
legend(strcat('n', num2str((1:nvars).')))


% Plot windowed variance of oscillator displacements.
Y_movvar = movvar(Y, [10, 0], 0, 2);
% 
% figure(2)
% plot(Y_movvar.')
% legend(strcat('n', num2str((1:nvars).')))

% % Compute potential energy of all springs feeding into an oscillator
% potentialEnergy = zeros(nvars, length(tSpan));
% for k = 1:length(tSpan)
%     p = [0; Y(:, k); 0];
%     r = repmat(p, 1, nvars + 2);
%     dist = r' - r;
%     sumPotentials = sum(K .* dist.^2 / 2, 2);
%     potentialEnergy(:, k) = sumPotentials(2:nvars+1);
% end
% 
% figure(2)
% plot(potentialEnergy.')
% legend(strcat('n', num2str((1:nvars).')))
% 
% 
% % Compute kinetic energy of every oscillator
% kineticEnergy = diff(Y, [], 2).^2 / 2;
% kineticEnergy = [zeros([nvars, 1]), kineticEnergy];
% 
% figure(3)
% plot(kineticEnergy.')
% legend(strcat('n', num2str((1:nvars).')))
% 
% % Compute entropy of every oscillator 
% windowLen = 200;
% entropy = zeros(nvars, length(tSpan) - windowLen + 1);
% for k = 1 : length(tSpan) - windowLen + 1
%     window = Y(:, k:k+windowLen-1);
%     entropy(:, k) = Entropy(window.');
% end
% 
% figure(4)
% plot(entropy.')
% legend(strcat('n', num2str((1:nvars).')))


% Compute distance from perturbed nodes using changepoints.
pad = 100;
obsIdx = 1:nvars;
%changepointThresh = 10;
%[pertOrders, pertResponseTimes] = ChngptPertOrders(Y, pertIdx, obsIdx, pertTimes, pad, changepointThresh)

% Compute distance from perturbed nodes using correlation coefficients.
corrThresh = 0.5;
observedY = Y;
%predPertValues = CorrPertValues(observedY, pertIdx, pertTimes, pertLength, pad, corrThresh)
[pertOrders, pertValues] = GetPertOrders(observedY, nvars, pertIdx, obsIdx, pertTimes, pertLength, 'corr', corrThresh, pad)

% Compute the true distance in the network.
obsIdx = 1:nvars;
truePertOrders = TruePertOrders(mat, pertIdx, obsIdx)
