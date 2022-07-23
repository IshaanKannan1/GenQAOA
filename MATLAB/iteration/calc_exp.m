function exp_val = calc_exp(p, params, D, sL, sR)
%     calc_exp computes the expected energy matrix element of the evolved
%     state for our four-driver QAOA. See report for details.
%     
%     Usage: expectation = calc_exp(p, params, D, sL, sR);
%     
%     Arguments:
%         p = QAOA depth
%         params = a px4 matrix specifying alpha, beta, gamma, delta in that
%         order, for each index 1 through p
%         D = vertex degree - 1
%         sL = left Pauli, integer 1-4, corresponding to {I, X, Y, Z} respectively
%         sR = same as sL, but for right Pauli

    A = params(:,1);
    A = [A;0;0;-flip(A)];

    z = gen_z(p);
    H_p = iterate(z, params, p, D);

    fbar = avg_f(z, params(:, 2:end), p);
    hl = comp_h_iterfun(z, sL, p);
    hr = comp_h_iterfun(z, sR, p);
    
    exp_val = 0;
    for iL = 1:length(z)
        exp_val = exp_val + H_p(iL) * fbar(iL) * hl(iL) ...
             * (exp(-1i * (z(iL, :).*z) * A).' * (H_p .* fbar .* hr));
    end
end


function z = gen_z(p)
    z = de2bi(0:2^(2*p+2)-1);
    z = z .* 2 - 1;
end
