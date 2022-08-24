function guesses = next_params_guess(prior, p_old, p_new)
    Ds = unique([prior.D]);
    params = cell(1, length(Ds));
    for i = 1:length(Ds)
        minpars = prior(i).finalpars([prior(i).energies] == ...
                min([prior(i).energies]));
        params(i) = {paramReduceQMC(cell2mat(minpars(1)), Ds(i))};
    end
    % two guesses: interpolation of a and d and continuation of b, c vs
    % interpolation of a, repetition of b, c, and d    

    x = 0:p_old-1;
    xnew = 0:(p_old-1)/(p_new-1):p_old-1;
    guesses = cell(1, length(Ds));
    
    % horizontally concat guesses for each D
    for i = 1:length(Ds)
        param = cell2mat(params(i));
        guess_a = transpose(interp1(x, param(:, 1).', xnew))
        guess_b = [param(:, 2); zeros(p_new-p_old, 1)]
        guess_c = [param(:, 3); zeros(p_new-p_old, 1)]
        guess_d = transpose(interp1(x, param(:, 4).', xnew))
        guess = [guess_a guess_b guess_c guess_d];

        guess_b = [param(:, 2); param(1:p_new-p_old, 2)];
        guess_c = [param(:, 3); param(1:p_new-p_old, 3)];
        guess_d = [param(:, 4); param(1:p_new-p_old, 4)];
        guess = [guess guess_a guess_b guess_c guess_d];
        guesses(i) = {guess};
    end
end
