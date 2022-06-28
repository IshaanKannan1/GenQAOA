function e = calc_exp(p, params, D, sL, sR)
%     calc_exp computes the expected energy matrix element of the evolved
%     state for our four-driver QAOA. See report for details.
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
    A = cat(1, A, -flip(A));
    z = gen_z(p);
    hps = iterate(z, params, p, D);
    fbar = transpose(avg_f(z, params(:, 2:end), p));
    hl = transpose(comp_h_iterfun(z, sL, p));
    hr = transpose(comp_h_iterfun(z, sR, p));
    e = 0;
    for i = 1:length(z)
        e = e + sum(hps(i) * fbar(i) * hl(i) * transpose(exp(-1i * (z(i).*z) * A)) .* ...
            hps .* fbar .* hr);
    end
end


function z = gen_z(p)
    z = de2bi(0:2^(2*p+2)-1);
    z = z .* 2 - 1;
end
