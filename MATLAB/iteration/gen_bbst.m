function [edges, A] = gen_bbst(depth)
% Genereates edge list and adjacency matrix of a graph described by two
% balanced binary trees joined at their roots, each of depth depth.
    A = zeros(2* (2^depth - 1));
    for j = 1:2^depth - 1 - 2^(depth - 1)
        A(j, 2*j:2*j+1) = 1;
        A(j + 2^depth - 1, 2^depth-1 +2*j:2^depth + 2*j) = 1;
    end

    A = A + transpose(A==1);
    A(1, 2^depth) = 1;
    A(2^depth, 1) = 1;
    [row, col] = find(triu(A));
    edges = [row, col];
end