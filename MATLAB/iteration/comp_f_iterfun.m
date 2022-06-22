function prod = comp_f_iterfun(n, z, params, p)
%     comp_f_iterfun efficiently computes each term in the iteration that 
%     factorizes over vertices, making use of some simplifications to avoid
%     repeated exponentiation. 
%     
%     Usage: out = comp_f_iterfun(randn, bitstring, params, p)
%     
%     Arguments:
%         n = a 3x1 vector on the unit sphere
%         z = a length 2p+2 bitstring representing our spin state basis.
%             0 represents a -1 in our algorithm.
%         params = a px3 matrix specifying parameters beta, gamma, delta
%                  in that order, for each depth from 1 to p
%         p = QAOA depth 

% compute product of exponential terms
    prod = 1;
    for i = 1:2*p+1
        R = z(i+1);
        L = z(i);
        beta_i = params(i, 1);
        gamma_i = params(i, 2);
        delta_i = params(i, 3);

        if i <= p
            if R == 1 && L == 1
                prod = prod * equal_z(beta_i, gamma_i, delta_i, n, 1);
            elseif R == -1 && L == -1
                prod = prod * equal_z(beta_i,gamma_i, delta_i, n, 0);
            elseif R == -1 && L == 1
                prod = prod * z1_pos_z2_neg(beta_i, gamma_i, delta_i, n);
            else
                prod = prod * z1_neg_z2_pos(beta_i, gamma_i, delta_i, n);
            end
        else
            if R == 1 && L == 1
                prod = prod * equal_z(beta_i, gamma_i, delta_i, n, 0);
            elseif R == -1 && L == -1
                prod = prod * equal_z(beta_i, gamma_i, delta_i, n, 1);
            elseif R == -1 && L == 1
                prod = prod * z1_pos_z2_neg_p(beta_i, gamma_i, delta_i, n);
            else
                prod = prod * z1_neg_z2_pos_p(beta_i, gamma_i, delta_i, n);
            end   
        end
    end
    z1 = z(1);
    zp = z(2*p+2);
    prod = prod * 1/4*(1+n(1)-1i*n(1)*zp+n(3)*zp + z1* ...
           (1i*n(2) + n(3) + zp - n(1)*zp));
end

function v = equal_z(b, g, d, n, pos)
    a = sqrt(b^2+2*b*d*n(1)+d^2*n(1)^2 + (g+d*n(2))^2 +d^2*n(3)^2);
    if pos == 1
        v = cos(a) + 1i*d*n(3)*sin(a)/a; 
    else
        v = cos(a) - 1i*d*n(3)*sin(a)/a;
    end
end

function v = z1_pos_z2_neg(b, g, d, n)
    a = sqrt(b^2+2*b*d*n(1)+d^2*n(1)^2 + (g+d*n(2))^2 +d^2*n(3)^2);
    v = (1i*b + g+ d*(1i*n(1) + n2)) * sin(a)/a;
end

function v = z1_neg_z2_pos(b, g, d, n)
    a = sqrt(b^2+2*b*d*n(1)+d^2*n(1)^2 + (g+d*n(2))^2 +d^2*n(3)^2);
    v = 1i*(b + 1i*g + d*(n(1) + 1i*n2)) * sin(a)/a;
end

function v = z1_pos_z2_neg_p(b, g, d, n)
    a = sqrt(b^2+2*b*d*n(1)+d^2*n(1)^2 + (g+d*n(2))^2 +d^2*n(3)^2);
    v = -1i*(b - 1i*(g + d*(1i*n(1) + n2))) * sin(a)/a;
end

function v = z1_neg_z2_pos_p(b, g, d, n)
    a = sqrt(b^2+2*b*d*n(1)+d^2*n(1)^2 + (g+d*n(2))^2 +d^2*n(3)^2);
    v = (-1i*b + g + d*(-1i*n(1) + n2)) * sin(a)/a;
end
