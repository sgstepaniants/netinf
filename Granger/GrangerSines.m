% Simulate data of two sine curves with (possibly) different amplitudes,
% phases, and frequencies. Apply Granger to this data and display the
% resulting causal connections predicted by MVGC.

clear all; close all; clc;

fn1 = @(x) 0.1 * sin(x);
fn2 = @(x) 0.1 * sin(x-pi/2);

measParam = 0.01;

tSpan = linspace(0, 4*pi, 200);
Y = [fn1(tSpan); fn2(tSpan)];
X = AddGaussianNoise(Y, measParam);

figure(1)
plot(X')

[connectivity, diagnostics] = DemoMVGC(X);

figure(2)
imagesc(connectivity)

figure(3)
load('A.mat')
p = size(A, 3);
ComputeVarFit(X, A, p)
