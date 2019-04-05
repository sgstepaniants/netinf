function dydt = odeNN(t, y, param)
%
% definition of ODE for NNCoupled model, to be passed to ode45
% 
% INPUTS:
%
% t
%       current time (not used)
%
% y
%       an [2n x 1] vector of the current state (y and y' stacked on top
%       of each other)
%
% param 
%       parameters that we need for the NNCoupled model
%       param{1} = number of nodes, n
%       param{2} = [n x 1] masses, m
%       param{3} = matrix of spring constants, K
%       param{4} = vector of damping coefficients, c
%       param{5} = boundary conditions
%
% OUTPUTS:
%
% dydt
%       an [2n x 1] vector of y' and y'' stacked on top of each other
    
    n = param{1};
    K = param{2};
    m = param{3};
    c = param{4};
    bc = param{5};
    f = param{6};
    
    p = y(1 : n); % p = [p1; p2; ...; p_n]
    p = p - uniffn(n, -0.5, 0.5, bc);
    v = y(n + 1 : 2 * n);

%     if strcmp(bc, 'circ')  % n x n springs for circular bc's
%         r = repmat(p, 1, n);
%         dist = r' - r;
%         force = sum((K + K') .* dist, 2);
%     elseif strcmp(bc, 'fixed')  % (n + 2) x (n + 2) springs for fixed bc's (wall anchored)
%         p = [0; p; 0];
%         r = repmat(p, 1, n + 2);
%         dist = r' - r;
%         force = sum(K .* dist, 2);
%         force = force(2 : (n + 1));
%     else  % n x n springs for free bc's (no springs on ends)
%         r = repmat(p, 1, n);
%         dist = r' - r;
%         force = sum(K .* dist, 2);
%     end
    
    p = [0; p; 0];
    r = repmat(p, 1, n + 2);
    dist = r' - r;
    force = sum(K .* dist, 2);
    force = force(2 : (n + 1));
    
    % calculate y''
    dydt = [v; (force - c .* v + f(t)) ./ m];
end
