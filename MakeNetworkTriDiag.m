function A = MakeNetworkTriDiag(nvars, stacks, corners)

% Make a random Erdos-Renyi network: for each pair of nodes,
% add a connection with probability p 
% (uniform distribution)
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
    stacks = 1;
    corners = true;
elseif nargin == 2
    corners = true;
end

block = ones(stacks);
N = nvars;

Dr = repmat(block, 1, N);
Dc = mat2cell(Dr, stacks, repmat(stacks, 1, N));
if ~corners
    Dc{N} = zeros(stacks);
end

D = blkdiag(Dc{:});
lowerBlockDiag = circshift(D, stacks);
upperBlockDiag = lowerBlockDiag.';

A = double(lowerBlockDiag | upperBlockDiag);
