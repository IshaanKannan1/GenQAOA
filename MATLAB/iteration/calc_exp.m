function e = calc_exp(p, params, D, sL, sR)
%     calc_exp computes the expected energy matrix element of the evolved
%     state for our four-driver QAOA. See our note for details.
%     
%     Usage: expectation = calc_exp(p, params, D, sL, sR);
%     
%     Arguments:
%         p = QAOA depth
%         params = a px4 matrix specifying alpha, beta, gamma, delta in that
%         order, for each index 1 through p
%         D = degree
%         sL = left Pauli, integer 1-4, corresponding to {I, X, Y, Z} respectively
%         sR = same as sL, but for right Pauli

    A = params(:,1);
    A(p+1) = 0;
    A(p+2) = 0;
    A = cat(1, A, flip(A));
    z = gen_z(p);
    e = 0;
    for i = 1:length(z)
        for j = 1:length(z)
            fl = avg_f(z(i), params(:, 2:end), p);
            fr = avg_f(z(j), params(:, 2:end), p);
            e = e + exp(-1i * dot(A, z(i).*z(j))) * iterate(z(i), fl, A, p, D) * ...
                    iterate(z(j), fr, A, p, D) * fl * fr * comp_h_iterfun(z(i), sL) * ...
                    comp_h_iterfun(z(j), sR);
        end
    end
end


function z = gen_z(p)
    a = dec2bin(0:2^p-1);
    z = (cell2mat(a) == '1');
    z = z .* 2 - 1;
end