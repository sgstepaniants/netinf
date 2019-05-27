% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

warning('off', 'stats:kmeans:EmptyCluster')

nvars = 3;
bc = 'fixed';

endtime = 200;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.1);

pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.2;
cfn = @(n) constfn(n, damping);

% Create connectivity matrix.
prob = 0.5;
spring = 0.1;
mat = MakeNetworkER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars + 2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = spring * K;


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
pertIdx = 1;
numPerts = length(pertIdx);
times = round(linspace(0, length(tSpan), numPerts+2));
pertTimes = times(2:end-1);
pertLength = 100;
pertForce = 10;
for k=1:numPerts
    forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
end

eps = 0.1;
minWaitTime = log(eps / min(sqrt(sum(((pertForce * inv(A)).^2))))) / max(amplitudes)
maxWaitTime = log(eps / max(sqrt(sum(((pertForce * inv(A)).^2))))) / max(amplitudes)

% Simulate harmonic oscillator movement.
ntrials = 1;
Y = GenerateHarmonicData(nvars, tSpan, 1, K, pfn, vfn, mfn, cfn, bc, forcingFunc);

% Plot oscillator trajectories.
figure(1)
plot(Y.'); hold on; plot((-inv(A) * forcingFunc).')
legend(strcat('n', num2str((1:nvars).')))


% Plot windowed variance of oscillator displacements.
% Y_movvar = movvar(Y, [10, 0], 0, 2);
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

leftPad = 500;
rightPad = pertLength;
obsIdx = true([1, nvars]);
%changepointThresh = 10;
%[pertOrders, pertResponseTimes] = ChngptPertOrders(Y, pertIdx, obsIdx, pertTimes, pad, changepointThresh)

% Compute distance from perturbed nodes using correlation coefficients.
corrThresh = 0.2;
observedY = noisefn(Y);
[pertOrders, pertValues] = GetPertOrders(observedY, pertIdx, obsIdx, pertTimes, leftPad, rightPad, 'corr', corrThresh)

% Compute the true distance in the network.
obsIdx = true([1, nvars]);
truePertOrders = TruePertOrders(mat, pertIdx, obsIdx)


%corrs = nan([nvars, 1000]);
%for t = 1 : 1000
%    corrMat = corr(observedY(:, (pertTimes - 500) : pertTimes + t).');
%    corrs(:, t) = corrMat(pertIdx, :);
%end
%figure(2);
%plot(corrs.');
