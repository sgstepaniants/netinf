clear all; close all; clc;
addpath('../DataScripts')
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

networkSizes = 1 : 20;
numSizes = length(networkSizes);

strengths = 1 : 25;
numStrengths = length(strengths);

numMats = 100;

endtime = 25;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);
prob = 0.5;

secondDerivThresh = 0.001;
w = gausswin(15);
w = w / sum(w);

charTimes = nan([numSizes, numStrengths, numMats]);
for j = 1 : numSizes
    nvars = networkSizes(j)
    for k = 1 : numStrengths
        strength = strengths(k)
        for m = 1 : numMats
            mat = MakeNetworkER(nvars, prob, true);
            forcingFunc = zeros(nvars, nobs);
            data = GenerateKuramotoData(mat, tSpan, 1, strength, pfn, wfn, forcingFunc);
            
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

save charTimes
meanCharTime = nanmean(charTimes, 3);

figure(1);
surf(meanCharTime)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca, 'TickLength', [0 0])


figure(2)
plot(1 ./ strengths(10:25), meanCharTime(10:end, 10:25), 'Linewidth', 2)
%hold on;
%plot(1 ./ strengths, 250 ./ strengths, 'Linewidth', 3)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca, 'TickLength', [0 0])
axis([0, 1.1, 0, 210])
