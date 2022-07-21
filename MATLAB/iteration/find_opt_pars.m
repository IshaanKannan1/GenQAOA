function [best_v, best_pars, best_D] = find_opt_pars(p, num_p, D)
    fname = join([fileparts(mfilename('fullpath')) strcat("/paramsearch/p:", num2str(p), ".mat")], "");
    % D = [1, 2, 3, 10, 100];
    data(num_p*length(D)) = struct();
    ctr = 1;
    for i = 1:length(D)
        for j = 1:num_p
            disp([D(i) j]);
            param0 = (2*rand(p, 4)-1) .* pi/2;
            [x, fval, exitflag] = opt_params(p, param0, D(i));
            data(ctr).D = D(i);
            data(ctr).energy = fval;
            data(ctr).exitflag = exitflag;
            data(ctr).initpars = param0;
            data(ctr).finalpars = mod(x, 2*pi);
            ctr = ctr + 1;
        end
    end
    save(fname, "data");
    [best_v, I] = min([data.energy]);
    best_v = -best_v;
    best_pars = data(I).finalpars;
    best_D = data(I).D;
end
