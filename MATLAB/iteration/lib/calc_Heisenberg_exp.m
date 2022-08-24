function [val, XX, YY, ZZ] = calc_Heisenberg_exp(p, params, D)
%calc_Heisenberg_exp computes the expected energy of the Heisenberg AF
%     Hamiltonian on a large-girth (D+1)-regular graph achieved by
%     the four-driver QAOA
%     
%     Usage: [val, XX, YY, ZZ] = calc_Heisenberg_exp(p, params, D);
%     
%     Arguments:
%         p = QAOA depth
%         params = a px4 matrix specifying alpha, beta, gamma, delta in 
%                  that order, for each index 1 through p
%         D = vertex degree - 1


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

    hX = comp_h_iterfun(zList, 2, p);
    hY = comp_h_iterfun(zList, 3, p);
    hZ = comp_h_iterfun(zList, 4, p);

    fH = H_p .* avg_f(zList, params(:, 2:end), p);
    fHX = fH.*hX; fHY = fH.*hY; fHZ = fH.*hZ;
    clear fH hX hY hZ

    batch_size = 2^p; % should divide 2^(2p+2)
    XX = 0;
    YY = 0;
    ZZ = 0;
    for ind = 1:batch_size:size(zList, 1)
        iL = ind-1 + (1:batch_size);

        temp = exp(-1i * (zList(iL, :).* A.') * zList.');
        XX = XX + fHX(iL).' * temp * fHX;
        YY = YY + fHY(iL).' * temp * fHY;
        ZZ = ZZ + fHZ(iL).' * temp * fHZ;
    end
    val = real(XX + YY + ZZ);
    
    % memory usage: (batch_size + 3) * 2^(2p+2)
end

