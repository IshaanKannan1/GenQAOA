function [vp, E] = d_inf_vp(p, params, varargin)

    A = params(:,1);
    A = [A; 0; 0; -flip(A)];
    
    zs = 2*de2bi(0:2^(2*p+2)-1) - 1;
    fs = avg_f(zs, params(:, 2:end), p);
    hs = cell(1, 4);
    hs{1} = comp_h_iterfun(zs, 1, p);
    hs{2} = comp_h_iterfun(zs, 2, p);
    hs{3} = comp_h_iterfun(zs, 3, p);
    hs{4} = comp_h_iterfun(zs, 4, p);

    G_curr = zs' * (fs .* zs); % does using sparse give a speedup?

    for i = 1:p-1 % go to p-1, because final vector needs p-1th matrix
        G_curr = get_next_G(G_curr, A, zs, fs, hs{1});
    end

    vp = cell(1, 3);
    m = exp(-0.5 * sum((zs * (A .* G_curr .* A')) .* zs, 2)); 
    for i = 1:3 
        F = fs.*hs{i + 1};
        Gp = zs' * (F .* m);
        vp{i} = 1i / 2 * sum(A .* (Gp.^2));
    end

    E = real((vp{1} + vp{2} + vp{3})/sqrt(varargin{1}));


    