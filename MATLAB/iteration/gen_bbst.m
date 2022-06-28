function [edges, A] = gen_bbst(depth)
% Tests iteration against actual QAOA result on a BBST
% Not great code but this will only run once
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