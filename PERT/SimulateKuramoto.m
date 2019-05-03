% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')

nvars = 5;

endtime = 10;
deltat = 0.01;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.1);

prob = 0.5;
mat = MakeNetworkER(nvars, prob, true); % adjacency matrix
K = 50; % connection strength

damping = 0;

pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - 1; % uniform [-1, 1]
cfn = @(n) constfn(n, damping);

% forcing function
pertForce = 50;
pertIdx = 1;
numPerts = length(pertIdx);
times = round(linspace(0, nobs, numPerts+2));
pertTimes = times(2:end-1);
pertLength = round(nobs/(3*(numPerts+1)));
forcingFunc = zeros([nvars, nobs]);
for k=1:numPerts
    forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
end

Y = GenerateKuramotoData(mat, tSpan, 1, K, pfn, wfn, cfn, forcingFunc);

figure(1)
plot(Y.')
legend(strcat('n', num2str((1:nvars).')))

% figure(2)
% plot(abs(diff(Y, [], 2)).')
% legend(strcat('n', num2str((1:nvars).')))
% 
% figure(3)
% plot(abs(diff(diff(Y, [], 2), [], 2)).')
% legend(strcat('n', num2str((1:nvars).')))

Y_movvar = movvar(Y, [10, 0], 0, 2);

figure(4)
plot(Y_movvar.')
legend(strcat('n', num2str((1:nvars).')))

% % Compute kinetic energy of every oscillator
% kineticEnergy = diff(Y, [], 2).^2 / 2;
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
% pad = 100;
% obsIdx = true([1, nvars]);
% changepointThresh = 10;
% [pertOrders, pertResponseTimes] = ChngptPertOrders(diff(Y, [], 2), pertIdx, obsIdx, pertTimes, pad, changepointThresh)

% Compute distance from perturbed nodes using correlation coefficients.
%Y_detrend = linearDetrend(Y, pertTimes);
%Y_detrend = detrend(Y.', 'linear', [pertTimes; pertTimes + pertLength]).';
%Y_detrend = Y_detrend - repmat(mean(Y_detrend(:, pertTimes-pad:pertTimes-1), 2), [1, length(tSpan)]);


windowData = Y(:, pertTimes:pertTimes+pertLength);
windowDataMovvar = movvar(windowData, [10, 0], 0, 2);
pertMeans = mean(windowDataMovvar, 2);
pertAve = pertMeans(pertIdx);
pertMeans(pertIdx) = [];
pertMeans = pertMeans / pertAve;

movvarWidth = nobs / 500;
obsIdx = true([1, nvars]);
predMeanPertOrders = MeanVarPertOrders(Y, pertIdx, obsIdx, pertTimes, pertLength, movvarWidth, 0)

% Compute the true distance in the network.
truePertOrders = TruePertOrders(mat, pertIdx, obsIdx)

figure(2); imagesc(truePertOrders); figure(3); imagesc(predMeanPertOrders)

% uu = normalize(truePertOrders);
% vv = normalize(-pertMeans);
% figure(3)
% plot(uu, 'g*'); hold on; plot(vv, 'r*')
