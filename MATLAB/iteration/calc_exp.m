function exp_val = calc_exp(p, params, D, sL, sR)
%     calc_exp computes the expected energy matrix element of the evolved
%     state for our four-driver QAOA. See report for details.
%     
%     Usage: exp_val = calc_exp(p, params, D, sL, sR);
%     
%     Arguments:
%         p = QAOA depth
%         params = a px4 matrix specifying alpha, beta, gamma, delta in 
%                  that order, for each index 1 through p
%         D = vertex degree - 1
%         sL = left Pauli, integer 1-4, corresponding to {I, X, Y, Z} respectively
%         sR = same as sL, but for right Pauli


    A = params(:,1);
    A = [A;0;0;-flip(A)];

    bList = logical(de2bi(0:2^(2*p)-1, 2*p));
    bList = [bList(:,1:p), false(2^(2*p),2), bList(:,p+1:end)];
    bList2 = bList; bList2(:,p+1) = true;
    bList3 = bList; bList3(:,p+2) = true;
    bList4 = bList2; bList4(:,p+2) = true;
    zList = 2*[bList;bList2;bList3;bList4] - 1;

    H_p = iterate(2*bList-1, params, p, D);
    H_p = [H_p; H_p; H_p; H_p];

    clear bList bList2 bList3 bList4

    fbar = avg_f(zList, params(:, 2:end), p);
    hl = comp_h_iterfun(zList, sL, p);
    hr = comp_h_iterfun(zList, sR, p);

    exp_val = 0;
    batch_size = 2^p; % should divide 2^(2p+2)

    for ind = 1:batch_size:length(zList)
        iL = ind-1 + (1:batch_size);

        exp_val = exp_val + sum( H_p(iL) .* fbar(iL) .* hl(iL) ...
             .* (exp(-1i * (zList(iL, :).* A.') * zList.' ) * (H_p .* fbar .* hr)) );
    end

end

