function H_p = iterate(z, params, p, D) 
% Computes all H_D(z) as outlined in our definition of the iteration. Returns
% a column of the p^th one. 
% 
% Usage: hp = iterate(z, fbar, A, p, D)
% 
% Arguments: 
%     z = array of bitstrings of length p of 1s and -1s.
%     param = px4 matrix of parameters, where the columns are in the order
%             of [alpha, beta, gamma, delta]
%     p = QAOA depth
%     D = degree - 1
%
% Returns:
%     H_D^{(p)} (z) for all the z, at the parameters specified


    A = params(:,1);
    A = [A;0;0;-flip(A)];

    xs = gen_all_x(p);
    fbar = avg_f(xs, params(:, 2:end), p);


    batch_size = 2^p; % any choice that divides 2^(2p+1) is ok

    H_prev = ones(2^(2*p+1), 1);
    H_next = H_prev;

    for i = 1:p-1
        
        % populate elements of H in batches
        for j = 1:batch_size:length(xs)
            Is = j-1 + (1:batch_size);
            
            % Note: fbar, H_prev are Nx1 vectors, where N=2^(2p+1)
            % xs = N x 2p+2, A = 2p+2 x 1, y = xs(Is, :) = L x 2p
            H_next(Is) = ((fbar.* H_prev).' * exp(-1i*(xs.*A.') * xs(Is, :).') ).^D;
        end
        H_prev = H_next;
    end

    H_p = zeros(size(z,1), 1);
    for ind = 1:batch_size:length(H_p)
        Is = ind-1 + (1:batch_size);
        H_p(Is) = ((fbar.* H_prev).' * exp(-1i*(xs.*A.') * z(Is, :).') ).^D;
    end
    Iremain = Is(end)+1:length(H_p);
    H_p(Iremain) =  ((fbar.* H_prev).' * exp(-1i*(xs.*A.') * z(Iremain, :).') ).^D;
    
% %   Alternative option to generate H_p, may not be memory-efficient:
%     H_p = transpose( ((fbar.* H_prev).' * exp(-1i*(xs.*A.') * z.') ).^D );

end
