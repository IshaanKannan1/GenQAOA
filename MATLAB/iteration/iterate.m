function hp = iterate(z, params, p, D) 
% Computes all H_D(z) as outlined in our definition of the iteration. Returns
% a column of the p^th one. 
% 
% Usage: hp = iterate(z, fbar, A, p, D)
% 
% Arguments: 
%     z = array of bitstrings of length p of 1s and -1s.
%     fbar = length p column of decimal values of fbar corre
%     A = 1x2p+2 array of alpha values, with 0s in the middle, as in note
%     p = QAOA depth
%     D = degree - 1

    A = params(:,1);
    A(p+1) = 0;
    A = cat(1, A, -flip(A));
    xs = gen_all_x(p);
    hds = zeros(p, 2^(2*p+1));
    hds(1, :) = 1;
    fbar = avg_f(xs, params(:, 2:end), p);
    for i = 2:p
         for j = 1:length(xs)
             hds(i, j) = sum(transpose(exp(-1i*(xs(j, :).*xs)*A)) .* hds(i-1, :) .* transpose(fbar)).^D;
             % hds(i, j) = (transpose(exp(-1i*(xs(j, :).*xs)*A)) .* hds(i-1, :) * fbar).^D;
%         for j = 1:2^p:length(xs)-2^p
%             Is = j-1 + (1:2^p);
%             % xs = N x 2p, A = 2p x 1, y = xs(Is, :) = L x 2p
%             hds(i, Is) = transpose(exp(-1i*(xs.*A.') * xs(Is, :).') .* hds(i-1, :) * fbar).^D;
        end
    end
    hp = zeros(1, size(z,1));
    for i = 1:length(hp)
        hp(i) = sum(transpose(exp(-1i*(z(i, :).*xs)*A)) .* hds(p, :) .* transpose(fbar)).^D;
    end
end
