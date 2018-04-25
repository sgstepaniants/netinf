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
%       param{3} = vector of spring constants, k
%       param{4} = vector of damping coefficients, c
%       param{5} = boundary conditions
%
% OUTPUTS:
%
% dydt
%       an [2n x 1] vector of y' and y'' stacked on top of each other

    global prev_pos;
    % dydt = vector of phase and its first derivative
    % y = phase
    
    % dydt(i, 2) = [k(i - 1) / m(i)] y(i - 1) - 2[k(i) / m(i)] y(i) + [k(i + 1) / m(i)] y(i + 1)
    n = param{1};
    m = param{2};
    k = param{3};
    c = param{4};
    bc = param{5};
    
    % reshape [y, y'] vector into two column matrix
    y = reshape(y, [n, 2]);
    % save previous positions of nodes
    prev_pos = [prev_pos, y(:, 1)];
    
    y0 = circshift(y(:, 1), [1, 0]);
    y1 = y(:, 1);
    y2 = circshift(y(:, 1), [-1, 0]);
    
    shift0 = y1 - y0;
    shift1 = y2 - y1;
    
    if strcmp(bc, 'circ')  % n springs for circular bc's
        if size(k, 1) == n
            k0 = circshift(k, [1, 0]);
            k1 = k;
        else
            error('spring constants for circular boundary conditions must have length n')
        end
        shift0(1) = shift0(1) + 1;
        shift1(n) = shift1(n) + 1;
    elseif strcmp(bc, 'fixed')  % n + 1 springs for fixed bc's (wall anchored)
        if size(k, 1) == n
            % if n spring constants given, make the spring constants of the
            % springs attached to the walls the same
            k0 = circshift(k, [1, 0]);
            k1 = k;
        elseif size(k, 1) == n + 1
            % if n + 1 spring constants given, set the prev spring constants
            % (k0) to the first n of these and set the next spring
            % constants (k1) to the last n of these.
            k0 = k(1 : n);
            k1 = k(2 : n + 1);
        else
            error('spring constants for fixed boundary conditions must have length n or n + 1')
        end
        shift0(1) = y(1, 1);
        shift1(n) = 1 - y(n, 1);
    else  % n - 1 springs for free bc's (no springs on ends)
        if size(k, 1) == n - 1
            % if n - 1 spring constants given, make outer spring constants
            % 0 (no springs)
            k0 = [0; k];
            k1 = [k; 0];
        elseif size(k, 1) == n
            % if n spring constants given, 0 out the last spring constant
            % because we only need n - 1 constants for free boundaries
            k(n) = 0;
            k0 = circshift(k, [1, 0]);
            k1 = k;
        else
            error('spring constants for free boundary conditions must have length n - 1 or n')
        end
        shift0(1) = 0;
        shift1(n) = 0;
    end
    
    % calculate y''
    dydt = zeros(n, 2);
    dydt(:, 1) = y(:, 2);
    dydt(:, 2) = (-k0 .* shift0 + k1 .* shift1 - c .* y(:, 2)) ./ m;
    
    % reshape [y', y''] two column matrix into vector
    dydt = reshape(dydt, [2 * n, 1]);
end
