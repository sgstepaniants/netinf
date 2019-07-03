clear all; close all; clc;

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

charTimes = nan([numSizes, numStrengths, numMats]);
for j = 1 : numSizes
    nvars = networkSizes(j)
    for k = 1 : numStrengths
        strength = strengths(k)
        for m = 1 : numMats
            mat = MakeNetworkER(nvars, prob, true);
            forcingFunc = zeros(nvars, nobs);

            data = GenerateKuramotoData(mat, tSpan, 1, strength, pfn, wfn, forcingFunc);
            secondDeriv = diff(diff(data, [], 2), 2);
            charTime = find(all(abs(secondDeriv) < 0.01, 1), 1, 'first');
            if isempty(charTime)
                charTime = NaN;
            end

            charTimes(j, k, m) = charTime;
        end
    end
end

meanCharTime = nanmean(charTimes, 3);
figure(1);
surf(meanCharTime)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca, 'TickLength', [0 0])


figure(2)
plot(1 ./ strengths, meanCharTime, 'Linewidth', 2)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca, 'TickLength', [0 0])
axis([0, 1.1, 0, 210])
