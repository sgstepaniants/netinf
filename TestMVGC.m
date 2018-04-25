load('UsualParams.mat')
reps = 10;

addpath('mvgc_v1.0');
startup;

% This test produces approximately 2 gigabiytes worth of experimental data.
% It runs GC on a nearest neighbor oscillator model with varying numbers of
% nodes, boundary conditions, initial conditions, and time splits.

nodes = 2 : 10;
bound_conds = {'free', 'fixed', 'circ'};
init_conds = {'rand', 'randstim', 'unifstim'};

% Specify fraction to split time into for voting.
fracs = [1, 0.5, 0.25, 0.125];

% Specify how much noise we want to add.
noisefn  = @(data) WhiteGaussianNoise(data, 0.12);
preprocfn = @(data) noisefn(data);

count = 1;
start_count = 73;
total_exp = size(nodes, 2) * size(bound_conds, 2) * size(init_conds, 2) * size(fracs, 2);

fileID = fopen('exp_log.txt','w');

for n = nodes
    for bound_cond = bound_conds
        if strcmp(bound_cond, 'circ')
            % Make the ground truth network with circular connections.
            mats = MakeNetworkTriDiag(n, 1, true);
        else
            % Make the ground truth network with single connections at
            % edge nodes.
            mats = MakeNetworkTriDiag(n, 1, false);
        end
        
        for init_cond = init_conds
            if strcmp(init_cond, 'rand')
                % Set the initial positions of the masses randomly.
                randpfn = @(n) rand(n, 1);
            elseif strcmp(init_cond, 'randstim')
                % Set the initial positions of the masses randomly.
                randpfn = @(n) rand(n, 1);
                % Set the intial velocity of the central node to 1.
                randvfn = @(n) [zeros(floor((n - 1) / 2), 1); 1; zeros(floor(n / 2), 1)];
            elseif strcmp(init_cond, 'unifstim')
                if strcmp(bound_cond, 'circ')
                    % Set the initial positions of the masses uniformly.
                    randpfn = @(n) [0 : n - 1]' / n;
                else
                    % Set the initial positions of the masses uniformly.
                    randpfn = @(n) [1 : n]' / (n + 1);
                end
                % Set the intial velocity of a central node to 1.
                randvfn = @(n) [zeros(floor((n - 1) / 2), 1); 1; zeros(floor(n / 2), 1)];
            end
            
            for frac = fracs
                % Get the experiment number.
                expnum = char(strcat(bound_cond, '_', init_cond, '_', ...
                                num2str(frac), '_', num2str(n)))
                
                % Calculate the times splits for this experiment.
                first = floor(frac * nobs);
                tsplits = first : first : nobs;
                
                % Run the base experiment.
                if (count >= start_count)
                    try
                        BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
                                        preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)
                    catch
                        fprintf(fileID, '***EXPERIMENT %d (%s) FAILED***', ...
                          count, expnum);
                    end
                end
                                
                % Count the number of experiments that have finished.
                fprintf(fileID, 'Experiment %d (%s) finished: %.1f%%\n', ...
                          count, expnum, count / total_exp * 100);
                count = count + 1;
            end
        end
    end
end

fclose(fileID);

exit;
