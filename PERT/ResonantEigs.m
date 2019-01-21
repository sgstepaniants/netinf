clear all; close all; clc;

%% For a fixed damping constant, plot statistics about getting resonant trials/eigenvalues

prob = 0.5;
numTrials = 1000;
numVars = 1000;

damping = 0.5;

freqAmplitudes = cell(numTrials, numVars);
for nvars = 1:numVars
    nvars
    for trial = 1:numTrials
        mat = MakeNetworkER(nvars, prob, true);
        K = MakeNetworkTriDiag(nvars + 2, false);
        K(2:nvars+1, 2:nvars+1) = mat;

        A = mat;
        A(1, 1) = A(1, 1) - 1;
        A(nvars, nvars) = A(nvars, nvars) - 1;
        A = A - diag(sum(mat, 2));
        lambdas = eig(A);
        normalfreqs = (-damping + sqrt(damping^2 + 4 * lambdas)) / 2;
        amplitudes = real(normalfreqs);
        
        freqAmplitudes{trial, nvars} = amplitudes;
    end
end


figure(1)
x = repelem(1:numVars, numTrials * (1 : numVars));
plot(x, vertcat(freqAmplitudes{:}).', 'b*')
xlabel('Number of Variables')
ylabel('Real Values of Frequencies')
title('Real Values of Frequencies for all Trials')


figure(2)
x = repelem(1:numVars, numTrials);
propResonantFreqs = zeros(1, numVars * numTrials);
for nvars = 1:numVars
    for trial = 1:numTrials
        propResonantFreqs(numTrials * (nvars - 1) + trial) = nnz(freqAmplitudes{trial, nvars} > 0.00001) / nvars;
    end
end
plot(x, propResonantFreqs, 'r*')
xlabel('Number of Variables')
ylabel('Proportion of Resonant Frequency Amplitudes')
title('Proportion of Resonant Frequency Amplitudes for all Trials')


figure(3)
x = 1:numVars;
aveResonantFreqs = zeros(1, numVars);
for nvars = 1:numVars
    for trial = 1:numTrials
        aveResonantFreqs(nvars) = aveResonantFreqs(nvars) + nnz(freqAmplitudes{trial, nvars} > 0.00001);
    end
    aveResonantFreqs(nvars) = aveResonantFreqs(nvars) / numTrials;
end
plot(x, aveResonantFreqs, 'r*')
xlabel('Number of Variables (N)')
ylabel('Average Number of Resonant Frequency Amplitudes for all Trials')
title('Average Number of Resonant Frequency Amplitudes in a Simulation with N Variables')


figure(4)
x = 1:numVars;
aveResonantTrials = zeros(1, numVars);
for nvars = 1:numVars
    for trial = 1:numTrials
        aveResonantTrials(nvars) = aveResonantTrials(nvars) + any(freqAmplitudes{trial, nvars} > 0.00001);
    end
    aveResonantTrials(nvars) = aveResonantTrials(nvars) / numTrials;
end
plot(x, aveResonantTrials, 'r*')
xlabel('Number of Variables (N)')
ylabel('Probability of Simulating Resonant Trial')
title('Probability of Simulating Resonant Trial in a Simulation with N Variables')


%% Plot the probability of getting a resonant trial over all numbers of variables and damping constants

prob = 0.5;
numTrials = 100;
numVars = 50;
dampingConsts = 0 : 0.01 : 1;

probResonantTrial = zeros(numVars, length(dampingConsts));
for nvars = 1:numVars
    nvars
    for i = 1:length(dampingConsts)
        damping = dampingConsts(i);
        for trial = 1:numTrials
            mat = MakeNetworkER(nvars, prob, true);
            K = MakeNetworkTriDiag(nvars + 2, false);
            K(2:nvars+1, 2:nvars+1) = mat;
            
            A = mat;
            A(1, 1) = A(1, 1) - 1;
            A(nvars, nvars) = A(nvars, nvars) - 1;
            A = A - diag(sum(mat, 2));
            lambdas = eig(A);
            normalfreqs = (-damping + sqrt(damping^2 + 4 * lambdas)) / 2;
            amplitudes = real(normalfreqs);
            
            probResonantTrial(nvars, i) = probResonantTrial(nvars, i) + any(amplitudes > 0.00001);
        end
        probResonantTrial(nvars, i) = probResonantTrial(nvars, i) / numTrials;
    end
end

figure(5)
imagesc(probResonantTrial.')
xlabel('Number of Variables')
ylabel('Damping Constant')
set(gca, 'XTick', 1:numVars, 'XTickLabel', 1:numVars)
set(gca, 'YTick', 1:length(dampingConsts), 'YTickLabel', dampingConsts)
