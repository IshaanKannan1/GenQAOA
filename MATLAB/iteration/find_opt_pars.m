function [best_v, best_pars] = find_opt_pars(p, num_p)
% Test D = 1, 2, 3, 100  
    fname = strcat("iteration/paramsearch/p:", num2str(p), ".mat");
    D = [1, 2, 3, 100];
    data = cell(length(D), num_p);
    best_v = 10000;
    best_pars = ones(1, 4);
    for i = 1:length(D)
        for j = 1:num_p
            param0 = (2*rand(p, 4) - 1) .* pi/2;
            [x, fval, exitflag] = opt_params(p, param0, D(i));
            data{i, j} = {param0, x, exitflag, fval};
            if fval < best_v
                best_v = fval;
                best_pars = x;                
            end
        end
    end
    save(fname, "data");
end