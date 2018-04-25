function Y = nncoupled(p, v, m, k, c, tSpan, bc, pert)
% calls ode45 to solve nearest neighbor coupled spring mass model
% 
% INPUTS:
%
% p
%       an [n x 1] vector of initial positions for the ODE
%
% v
%       an [n x 1] vector of initial velocities for the ODE
%
% adj
%       an [n x n] adjacency matrix for the network
%
% m 
%       an [n x 1] vector of masses, one for each node
%       in the network
%
% k 
%       an n - 1, n, or n + 1 length vector of spring constants,
%       one for each edge in the network
%
% c 
%       an n - 1, n, or n + 1 length vector of damping coefficients,
%       one for each edge in the network
%
% tSpan
%       a [1 x m] vector of times for which you want the solution of the
%       ODE returned
%
% bc
%       a string that specifies the boundary conditions of the nearest
%       neighbor spring model ('circ', 'fixed', 'free')
%
% pert  a number that signifies how much nodes should be perturbed when
%       their oscillation amplitudes fall below a certain threshold
%
% OUTPUTS:
%
% Y
%       an [n x m] matrix of the solution of the ODE for all n 
%       nodes and m time points
    global prev_pos;
    global thresh;
    global past;

    n = size(p, 1);
    param{1} = n;     % number of nodes
    param{2} = m;     % masses
    param{3} = k;     % spring constants
    param{4} = c;     % damping coefficients
    param{5} = bc;    % boundary conditions

    tstart = tSpan(1);
    tfinal = tSpan(end);
    y0 = [p; v];
    options = odeset('Events', @events);
    
    % Solve until the first terminal event
    sol = ode45(@(t, y) odeNN(t, y, param), [tstart tfinal], y0, options);
    
    while sol.x(end) < tfinal
        %fprintf('"%i is where it actually terminates"\n', sol.x(end))
        % Find which nodes have amplitudes below a certain threshold
        amp = range(prev_pos(:, end - past + 1 : end), 2) / 2;
        
        % Set the new initial conditions, with perturbed nodes
        p0 = sol.y(1:n, end);
        v0 = sol.y(n+1:2*n, end) - pert * randn(n, 1) .* (amp < thresh) ./ m;
        y0 = [p0; v0];
        
        % Accumulate output
        sol = odextend(sol, @(t, y) odeNN(t, y, param), tfinal, y0, options);
    end
    
    Y = deval(sol, tSpan);
    
    % take only the top portion of the matrix which holds the positions
    % of each one of the n masses over time
    Y = Y(1:n, :, :);

%---------------------------------------------------------------------------

function [value, isterminal, direction] = events(t, y)
% If some of the nodes have amplitudes smaller than a certain threshold,
% stop the integration
    global prev_pos;
    global thresh;
    global past;
    
    % Get the number of nodes in the system
    n = size(y, 1) / 2;
    
    value = ones(n, 1);
    if size(prev_pos, 2) >= past
        amp = range(prev_pos(:, end - past + 1 : end), 2) / 2;
        value = double(amp >= thresh);  % detect amplitude < thresh
        
        % mess with the recorded previous node positions to ensure that
        % the nodes are not kicked on two consecutive iterations
        prev_pos(:, end - past + 20) = prev_pos(:, end) - 4 * amp;
    end
    if any(value)
        %fprintf('%i is where it should terminate\n', t)
    end
    isterminal = ones(n, 1);  % stop the integration
    direction = zeros(n, 1);  % no direction
