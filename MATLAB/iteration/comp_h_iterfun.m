function v = comp_h_iterfun(z, sigma, p)
% Computes small h function containing middle term, based on desired Pauli
% or identity. Does not use any matrix multiplication for efficiency.
% 
% Arguments:
%     z = matrix of bitstrings of 1s and -1s as array of length 2p+2
%     sigma = {1, 2, 3, 4} corresponding to {I, X, Y, Z} respectively.
    z1 = z(:, p+1);
    z2 = z(:, p+2);
    switch sigma
        case 1
            v = 1/2*(1+z1.*z2);
        case 2
            v = 1/2*(1-z1.*z2);
        case 3
            v = -1i/2*(z1-z2);
        case 4
            v = (z1 + z2)/2;
    end
end