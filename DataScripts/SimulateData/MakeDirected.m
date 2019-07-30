function directedMat = MakeDirected(mat, keepConn)
    if nargin < 2
        keepConn = 0;
    end
    
    n = size(mat, 1);
    deleteEdges = randi([0, 1], n);
    if keepConn
        deleteEdges = tril(deleteEdges, -1) + ~tril(deleteEdges, -1).';
    end
    directedMat = mat .* deleteEdges;
end
