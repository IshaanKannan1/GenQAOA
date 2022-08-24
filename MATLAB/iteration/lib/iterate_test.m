%% This file tests the subroutines in iterate(), used to generate H_D^(m)
%% Setting up input of iterate()
p = 6;

params = ones(p,4);
D = 2;
z = 2*(de2bi(0:2^(2*p+2)-1) - 1);

%%

A = params(:,1);
A = [A;0;0;-flip(A)];

xs = gen_all_x(p);

fbar = avg_f(xs, params(:, 2:end), p);


%% old version (before 7/22/2022)
tic;
    hds = ones(p+1, 2^(2*p+1));

    for i = 2:p
         for j = 1:length(xs)
             hds(i, j) = sum(transpose(exp(-1i*(xs(j, :).*xs)*A)) .* hds(i-1, :) .* transpose(fbar)).^D;
        end
    end
    H_p_old = zeros(size(z,1), 1);
    for i = 1:length(H_p_old)
        H_p_old(i) = sum(transpose(exp(-1i*(z(i, :).*xs)*A)) .* hds(p, :) .* transpose(fbar)).^D;
    end

toc;

%% new version (after 7/22/2022)
tic;

    batch_size = 2^(p); % any choice that divides 2^(2p+1) is ok

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
    
%     H_p = transpose( ((fbar.* H_prev).' * exp(-1i*(xs.*A.') * z.') ).^D );

toc;