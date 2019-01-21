function Y = kuramoto(p,adj,w,c,tSpan,K, forcingFunc)
% calls ode45 to solve kuramoto model
% 
% INPUTS:
%
% p
%       an [n x 1] vector of initial conditions for the ODE
%
% adj
%       an [n x n] adjacency matrix for the network
% 
% w 
%       an [n x 1] vector of natural frequencies, one for each node
%       in the network
%
% c 
%       an [n x 1] vector of damping coefficients, one for each node
%       in the network
%
% tSpan
%       a [1 x m] vector of times for which you want the solution of the
%       ODE returned
%
% K
%       a scalar value for connection strength of the network 
%
% OUTPUTS:
%
% Y
%       an [n x m] matrix of the solution of the ODE for all n 
%       nodes and m time points
    
    n = size(p,1);
    
    param{1} = n; % number of nodes
    param{2} = adj; % adjacency matrix
    param{3} = w; % natural frequencies
    param{4} = K; % connection strength
    param{5} = c; % damping coefficients
    
    if nargin < 7
        forcingFunc = zeros(n, length(tSpan));
    end
    param{6} = @(t) forcingFunc(:, findClosestTime(t, tSpan));
    
    % Solve until the first terminal event
    SOL = ode45(@(t, y) odeKur(t, y, param), tSpan, p);
    Y = deval(SOL, tSpan);
end

function idx = findClosestTime(t, tSpan)
    [~, idx] = min(abs(tSpan-t));
end
