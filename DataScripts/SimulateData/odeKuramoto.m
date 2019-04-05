function dy = odeKuramoto(t,y,param)
%
% definition of ODE for Kuramoto model, to be passed to ode45
% 
% INPUTS:
%
% t
%       current time (not used) 
%
% y
%       an [2n x 1] vector of the current state
% 
% param 
%       parameters that we need for the Kuramoto model
%       param{1} = number of nodes, n
%       param{2} = [n x n] adjacency matrix, A
%       param{3} = [n x 1] natural frequencies, w
%       param{4} = scalar connection strength, K
%       param{5} = [n x 1] damping coefficients, c
%
% OUTPUTS:
%
% dy
%       an [n x 1] vector of y'


    % dy = angular velocity
    % y  = phase
    n = param{1};
    A = param{2};
    w = param{3};
    K = param{4};
    c = param{5};
    f = param{6};
    
    % dp = w + sum over j (k * sin(y(j) - y(i)))
    r = repmat(y,1,n);
    
    % adj matrix: Aij non-zero if node j (colm) infl. node i (row)
    % sum over the columns (all the influences) 
    % ijth entry of r' - r should be jth node - ith node
    % example: 5,1 entry is 1st node - 5th node
    dy = w + (K/n)*sum(A .* sin(r'-r),2) - c.*y + f(t);
end
