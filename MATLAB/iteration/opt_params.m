function [x, fval, exitflag] = opt_params(p, params, D)  
    optfun = @(param0) 1/2 * (1 - calc_exp(p, param0, D, 2, 2) - ...
                             calc_exp(p, param0, D, 3, 3) - calc_exp(p, param0, D, 4, 4));
    options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
        'TolX',1e-5,'TolFun',1e-5, 'Algorithm', 'quasi-newton','PlotFcns',...
        'MaxFunEvals', Inf, 'MaxIter', Inf);
    [x, fval, exitflag] = fminunc(optfun, params, options);
%     val = fval;
%     paramsf = x;
%     exit = exitflag;
end