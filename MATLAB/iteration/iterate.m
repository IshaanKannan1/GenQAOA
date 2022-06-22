function hp = iterate(z, fbar, A, p, D) 
% Computes all H_D(z) as outlined in our definition of the iteration. Returns
% the p^th one. 
% 
% Usage: hp = iterate (z, fbar, A, p, D)
% 
% Arguments: 
%     z = any bitstring of length p of 1s and -1s provided in array form.
%     fbar = decimal value of fbar
%     A = 1x2p+2 array of alpha values, with 0s in the middle, as in note
%     p = QAOA depth
%     D = degree

    xs = gen_all_x(p);
    hds = zeros(p, 1);
    hds(1) = 1;
    for i = 2:p
        for j = 1:length(xs)
            hds(i) = hds(i) + (exp(-1i*dot(A, xs(j).* z))) * hds(i-1) * fbar;
        end
        hds(i)=hds(i)^D;
    end
    hp = hds(p);
end
