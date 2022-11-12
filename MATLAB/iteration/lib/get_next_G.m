function G_next = get_next_G(G, A, zs, fs, hs)
    m = exp(-0.5 * sum((zs * (A .* G .* A')) .* zs, 2));  
    F = (fs.*hs);
    G_next = zs' * ((F .* m) .* zs);