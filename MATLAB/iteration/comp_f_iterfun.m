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
    prod = ones(size(z, 1), 1);
    j = 1;
    cos_z = cos(params(:, 1)) .* cos(params(:, 2));
    m = normr([cos(params(:, 2)).*sin(params(:, 1)), sin(params(:, 1)) .* sin(params(:, 2)), ...
         cos(params(:, 1)).*sin(params(:, 2))]);
    sin_z = sqrt(1-cos_z.^2);
    cos_t = cos(params(:, 3)) .* cos_z - sin_z.*sin(params(:, 3)) .* (m * n.')
    sin_t = sqrt(1-cos_t.^2)
    r = normr(sin_t .* cos(params(:, 3)) .* m + cos_t .* sin(params(:, 3)) .* n - sin_t .* sin(params(:, 3)) ...
        .* cross(m, repmat(n, p, 1)))
    while j <= 2*p+1
        R = z(:, j+1);
        L = z(:, j);
        if j <= p
            prod = prod .* ((cos_t(j) + 1i * sin_t(j) * r(j, 3)) * (R+1).*(L+1)/4 ...
                          + (cos_t(j) - 1i * sin_t(j) * r(j, 3))* (1-R).*(1-L)/4 ...
                          + ((r(j, 2) + 1i * r(j, 1)) * sin_t(j)) * (L+1).*(1-R)/4 ...
                          + ((-r(j, 2) + 1i * r(j, 1)) * sin_t(j)) * (1-L).*(R+1)/4);
        else
            k = 2*p+2 - j;
            prod = prod .* ((cos_t(k) - 1i * sin_t(k) * r(k, 3)) * (R+1).*(L+1)/4 ...
                          + (cos_t(k) + 1i * sin_t(k) * r(k, 3))* (1-R).*(1-L)/4 ...
                          - ((r(k, 2) + 1i * r(k, 1)) * sin_t(k)) * (L+1).*(1-R)/4 ...
                          + ((r(k, 2) - 1i * r(k, 1)) * sin_t(k)) * (1-L).*(R+1)/4);
        end
% 
%         if j <= p
%             beta_i = params(j, 1);
%             gamma_i = params(j, 2);
%             delta_i = params(j, 3);
%             prod = prod .* (equal_z(beta_i, gamma_i, delta_i, n, 1) * (R+1).*(L+1)/4 ...
%                          + equal_z(beta_i,gamma_i, delta_i, n, 0) * (1-R).*(1-L)/4 ...
%                          + z1_pos_z2_neg(beta_i, gamma_i, delta_i, n) * (L+1).*(1-R)/4 ...
%                          + z1_neg_z2_pos(beta_i, gamma_i, delta_i, n) * (1-L).*(R+1)/4);
%         else
%             beta_i = params(2*p+2 - j,1);
%             gamma_i = params(2*p+2 - j, 2);
%             delta_i = params(2*p+2 - j, 3);
%             prod = prod .* (equal_z(-beta_i, -gamma_i, -delta_i, n, 1) * (R+1).*(L+1)/4 ...
%                          + equal_z(-beta_i, -gamma_i, -delta_i, n, 0) * (1-R).*(1-L)/4 ...
%                          + z1_pos_z2_neg(-beta_i, -gamma_i, -delta_i, n) * (L+1).*(1-R)/4 ...
%                          + z1_neg_z2_pos(-beta_i, -gamma_i, -delta_i, n) * (1-L).*(R+1)/4);
%         end
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

% 
% function v = equal_z(b, g, d, n, pos)
%     a = sqrt(b^2 + 2*b*d*n(1) + d^2*(n(1)^2 + n(2)^2) + (g + d*n(3))^2);
%     if pos == 1
%         v = exp(-1i*a) * (exp(2*1i*a)*(g + d*n(3) + a) - g - d*n(3) + a) / (2*a); 
%     else
%         v = exp(-1i*a) * (-exp(2*1i*a)*(g + d*n(3) - a) + g + d*n(3) + a) / (2*a);
%     end
% end
% 
% function v = z1_pos_z2_neg(b, g, d, n)
%     a = sqrt(b^2 + 2*b*d*n(1) + d^2*(n(1)^2 + n(2)^2) + (g + d*n(3))^2);
%     v = (1i*(b + d*n(1)) + d*n(2)) * sin(a)/a;
% end
% 
% function v = z1_neg_z2_pos(b, g, d, n)
%     a = sqrt(b^2+2*b*d*n(1)+d^2*n(1)^2 + (g+d*n(2))^2 +d^2*n(3)^2);
%     v = 1i*(b + d*(n(1) + 1i*n(2))) * sin(a)/a;
% end
