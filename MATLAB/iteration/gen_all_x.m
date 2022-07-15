function strs = gen_all_x(p)
% Helper function to generate array of all bitstrings of length 2p+2 where 
% two centermost bits are equal (all 1s or -1s). Returns 2^(2p+1)x2p+2
% matrix of 1s and -1s. 
    a = de2bi(0:2^p-1);
    l = [[0, 0]; [1, 1]];
    strs = zeros(2^(2*p+1),2*p+2);
    ctr = 1;
    for j = 1:2
        for k = 1:2^p
            for i = 1:2^p
                strs(ctr, :) = [a(k,:), l(j,:), a(i,:)];
                ctr = ctr + 1;
            end
        end
    end
    strs = strs .* 2 - 1;
end