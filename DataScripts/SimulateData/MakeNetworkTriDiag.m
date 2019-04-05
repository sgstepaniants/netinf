function A = MakeNetworkTriDiag(nvars, corners, stacks)

% Make a tridiagonal network
%
% INPUTS:
%
% nvars
%       scalar for number of nodes in network
%
% corners
%       boolean whether the matrix should have values in the corners
%       for a system with circular boundary conditions
%
% OUTPUTS:
%
% A
%       adjacency matrix for the network created
%

if nargin == 1
    corners = true;
    stacks = 1;
elseif nargin == 2
    stacks = 1;
end

block = ones(stacks);

Dr = repmat(block, 1, nvars);
Dc = mat2cell(Dr, stacks, repmat(stacks, 1, nvars));
if ~corners
    Dc{nvars} = zeros(stacks);
end

D = blkdiag(Dc{:});
lowerBlockDiag = circshift(D, stacks);
upperBlockDiag = lowerBlockDiag.';

A = double(lowerBlockDiag | upperBlockDiag);
