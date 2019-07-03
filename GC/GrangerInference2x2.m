% Show Granger works on a 2x2 Network

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

nvars = 2;

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

numTrialsList = 5;
numTrialsListLength = length(numTrialsList);

numMats = 100;

trueMats = zeros(nvars, nvars, numTrialsListLength, numMats);
predMats = zeros(nvars, nvars, numTrialsListLength, numMats);
accLog = zeros(numTrialsListLength, numMats);
for j = 1:numTrialsListLength
    numTrials = numTrialsList(j)
    
    for m = 1:numMats
        mat = MakeNetworkER(nvars, 0.5, true);
        trueMats(:, :, j, m) = mat;
        
        K = zeros(nvars + 2);
        K(1, 2) = 1;
        K(2, 1) = 1;
        K(end - 1, end) = 1;
        K(end, end - 1) = 1;
        K(2:nvars+1, 2:nvars+1) = mat;
        
        data = GenerateHarmonicData(nvars, tSpan, numTrials, K, ...
            pfn, vfn, mfn, cfn, bc, forcingFunc);
        
        predMat = DemoMVGC(preprocfn(noisefn(data)), 1);
        predMats(:, :, j, m) = predMat;
        accLog(j, m) = nnz((predMat == mat) .* ~eye(nvars)) / (nvars^2-nvars);
    end
end

acc = mean(accLog, 2);
plot(numTrialsList, acc)
xlabel('Number of Trials fed to GC')
ylabel('Accuracy')
