% Show Granger works on a 2x2 Network

clear all; close all; clc;

nvars = 20;

mfn = @(n) constfn(n, 1);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

preprocfn = @(data) standardize(data);

damping = 0.2;
cfn = @(n) constfn(n, damping);

deltat = 0.1;
endtime = 25;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

forcingFunc = zeros([nvars, nobs]);

bc = 'fixed';

numTrialsList = 1:100;

mat = MakeNetworkSymmER(nvars, 0.5, true);
K = zeros(nvars + 2);
K(1, 2) = 1;
K(2, 1) = 1;
K(end - 1, end) = 1;
K(end, end - 1) = 1;
K(2:end-1, 2:end-1) = mat;

data = GenerateHarmonicData(nvars, tSpan, max(numTrialsList), K, pfn, vfn, mfn, cfn, bc, forcingFunc);

numMats = 100;
predMats = zeros(nvars, nvars, numMats, length(numTrialsList));
accLog = zeros(length(numTrialsList), numMats);
for j = 1:length(numTrialsList)
    for k = 1:numMats
        numTrials = numTrialsList(j);
        predMat = DemoMVGC(preprocfn(noisefn(data(:, :, 1:numTrials))), 1);
        predMats(:, :, j, k) = predMat;
        accLog(j, k) = nnz((predMat == mat) .* ~eye(nvars)) / (nvars^2-nvars);
    end
end

acc = mean(accLog, 2);
plot(acc)
xlabel('Number of Trials fed to GC')
ylabel('Accuracy')
