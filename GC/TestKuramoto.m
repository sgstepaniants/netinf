% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;

nvars = 2;

endtime = 25;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);

prob = 0.5;
%mat = ones(nvars); % adjacency matrix
%mat(1 : nvars + 1 : nvars^2) = 0;
mat = MakeNetworkER(nvars, prob, true);
K = 10; % connection strength

pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

forcingFunc = zeros(nvars, nobs);
Y = GenerateKuramotoData(mat, tSpan, 1, K, pfn, wfn, forcingFunc);
X = noisefn(Y);
hax = axes;
plot(Y.')
%plot3(X(1, :), X(2, :), X(3, :), 'ro')

secondDerivThresh = 0.001;
w = gausswin(15);
w = w / sum(w);
smoothedData = filter(w, 1, Y, [], 2);
secondDeriv = diff(smoothedData, 2, 2);
trajAligned = all(abs(secondDeriv) < secondDerivThresh, 1);
charTime = find(diff(trajAligned), 1, 'last');
if trajAligned(end) == 0
    charTime = NaN;
elseif all(trajAligned) == 1
    charTime = 1;
end

hold on;
line([charTime, charTime], get(hax,'YLim'), 'Linewidth', 3, 'Color', 'k')

% for t = 1 : size(tSpan, 2)
%     for i = 1 : nvars
%         pos = 2 * pi * Y(i, t, 1);
%         r = 0.1;
%         rectangle('Position', [cos(pos) - r, sin(pos) - r, 2 * r, 2 * r], ...
%             'FaceColor',[i / nvars, 0.2, 0.2], 'Curvature', [1, 1])
%         axis([-2 2 -2 2])
%     end
%     pause(0.001)
%     clf
% end
