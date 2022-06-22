function strs = gen_all_x(p)
%     Generate array of all bitstrings of length 2p+2 where two centermost bits
%     are equal (all 1s or -1s). Returns 2^(2p+1)xp matrix of 1s and -1s.
   
    a = dec2bin(0:2^p-1);
    l = ['00'; '11'];
    strs = cell(2^(2*p+1),1);
    strs(:) = {''};
    ctr = 1;
    % is there a more efficient way?
    for j = 1:2
        for k = 1:2^p
            for i = 1:2^p
                strs(ctr) = {strcat(a(k,:), l(j,:), a(i,:))};
                ctr = ctr + 1;
            end
        end
    end
    strs = (cell2mat(strs) == '1');
    strs = strs .* 2 - 1;
end