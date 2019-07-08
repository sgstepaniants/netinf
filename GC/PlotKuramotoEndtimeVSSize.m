clear all; close all; clc;
addpath('../DataScripts')
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PlotEndtimeSize';
% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);

networkSizes = 1 : 20;
numSizes = length(networkSizes);

strengths = 1 : 25;
numStrengths = length(strengths);

numMats = 100;

endtime = 25;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);
prob = 0.5;

secondDerivThresh = 0.001;
w = gausswin(15);
w = w / sum(w);

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));


dataLog = cell(numSizes, numStrengths, numMats);
charTimes = nan([numSizes, numStrengths, numMats]);
for j = 1 : numSizes
    nvars = networkSizes(j)
    for k = 1 : numStrengths
        strength = strengths(k)
        for m = 1 : numMats
            mat = MakeNetworkER(nvars, prob, true);
            forcingFunc = zeros(nvars, nobs);
            data = GenerateKuramotoData(mat, tSpan, 1, strength, pfn, wfn, forcingFunc);
            dataLog{j, k, m} = data;
            
            smoothedData = filter(w, 1, data, [], 2);
            secondDeriv = diff(smoothedData, 2, 2);
            trajAligned = all(abs(secondDeriv) < secondDerivThresh, 1);
            charTime = find(diff(trajAligned), 1, 'last');
            if trajAligned(end) == 0
                charTime = NaN;
            elseif all(trajAligned) == 1
                charTime = 1;
            end

            charTimes(j, k, m) = charTime;
        end
    end
end

save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/charTimes.mat', expPath), 'charTimes');

meanCharTime = nanmean(charTimes, 3);
figure(1);
surf(meanCharTime)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca, 'TickLength', [0 0])


figure(2)
strengthRange = 1:25;
sizeRange = 1:20;
plot(1 ./ strengths(strengthRange), meanCharTime(sizeRange, strengthRange), 'Linewidth', 2)
%p = polyfit(1 ./ strengths(strengthRange), meanCharTime(sizeRange, strengthRange), 1);
hold on;
%plot(1 ./ strengths, p(0) ./ strengths + p(1), 'Linewidth', 3)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca, 'TickLength', [0 0])
axis([0, 1, 0, 250])
