% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;

nvars = 3;

endtime = 15;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);

prob = 0.5;
mat = ones(nvars); % adjacency matrix
mat(1 : nvars + 1 : nvars^2) = 0;
%mat = MakeNetworkER(nvars, prob, true);
K = 0.5; % connection strength

pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]

forcingFunc = zeros(nvars, nobs);
Y = GenerateKuramotoData(mat, tSpan, 1, K, pfn, wfn, forcingFunc);
X = noisefn(Y);
%plot(Y.')

plot3(X(1, :), X(2, :), X(3, :), 'ro')

secondDeriv = diff(diff(Y, [], 2), 2);
charTime = find(all(abs(secondDeriv) < 0.01, 1), 1, 'first');
if isempty(charTime)
    charTime = Inf;
end

charTime

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
