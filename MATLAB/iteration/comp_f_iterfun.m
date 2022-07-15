function prod = comp_f_iterfun(n, z, params, p)
%     comp_f_iterfun efficiently computes each term in the iteration that 
%     factorizes over vertices, making use of some simplifications to avoid
%     repeated exponentiation. 
%     
%     Usage: out = comp_f_iterfun(randn, bitstring, params, p)
%     
%     Arguments:
%         n = a 3x1 vector on the unit sphere
%         z = a matrix of any number of length 2p+2 bitstring representing 
%             our spin state basis.
%         params = a px3 matrix specifying parameters beta, gamma, delta
%                  in that order, for each depth from 1 to p
%         p = QAOA depth 
%
%         prod = f value for each bitstring, column of values

    betas = params(:, 1);
    gammas = params(:,2);
    deltas = params(:,3);
    cos_z = cos(betas) .* cos(gammas);

    sin_z = sqrt(1-cos_z.^2);
    m = [cos(gammas).*sin(betas), sin(betas) .* sin(gammas), ...
         cos(betas).*sin(gammas)] ./sin_z;
 
    cos_t = cos(deltas) .* cos_z - sin_z.*sin(deltas) .* (m * n.');
    sin_t = sqrt(1-cos_t.^2);

    r = (sin_z .* cos(deltas) .* m + cos_z .* sin(deltas) .* n - sin_z .* sin(deltas) ...
        .* cross(m.', repmat(n, p, 1).').' ) ./sin_t;

    prod = ones(size(z, 1), 1);

    j = 1;
    while j <= 2*p+1
        if j <= p
            R = z(:, j+1);
            L = z(:, j);
            prod = prod .* ( (cos_t(j) + 1i * sin_t(j) * r(j, 3)) * (1+L).*(1+R)/4 ...
                           + (cos_t(j) - 1i * sin_t(j) * r(j, 3)) * (1-L).*(1-R)/4 ...
                           + (1i * (r(j, 1) + 1i * r(j, 2)) * sin_t(j)) * (1+L).*(1-R)/4 ...
                           + (1i * (r(j, 1) - 1i * r(j, 2)) * sin_t(j)) * (1-L).*(1+R)/4);
        else
            k = 2*p+2 - j;
            R = z(:, j);
            L = z(:, j+1);
            
            prod = prod .* ...
                        conj((cos_t(k) + 1i * sin_t(k) * r(k, 3)) * (1+L).*(1+R)/4 ...
                           + (cos_t(k) - 1i * sin_t(k) * r(k, 3)) * (1-L).*(1-R)/4 ...
                           + (1i * (r(k, 1) + 1i * r(k, 2)) * sin_t(k)) * (1+L).*(1-R)/4 ...
                           + (1i * (r(k, 1) - 1i * r(k, 2)) * sin_t(k)) * (1-L).*(1+R)/4);
        end
    
        if j == p
            j = j + 2;
        else
            j = j + 1;
        end
    end

    z1 = z(:, 1);
    zp = z(:, 2*p+2);
    prod = prod .* 1/4.*(1 + n(1) + n(3).*z1 + 1i*n(2).*(z1-zp) - n(1).*zp.*z1 + (n(3) + z1).*zp);
end
