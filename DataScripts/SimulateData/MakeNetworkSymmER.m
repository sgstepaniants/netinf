function A = MakeNetworkSymmER(nvars, p, zeroDiag)

% Make a random symmetric Erdos-Renyi network: for each pair of nodes,
% add a connection with probability p 
% (uniform distribution)
%
% INPUTS:
%
% nvars
%       scalar for number of nodes in network
%
% p
%       probability than an edge is included in the network
%
% OUTPUTS:
%
% A
%       adjacency matrix for the network created
%

if nargin == 2
    zeroDiag = true;
end

er = MakeNetworkER(nvars, p, zeroDiag);
lowerTri = tril(er, -1);
A = lowerTri + lowerTri' + diag(diag(er));
