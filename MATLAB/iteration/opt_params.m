function [x, fval, exitflag, funcCount] = opt_params(p, params, D)  
%opt_params finds optimal 4-driver QAOA parameters for minimizing
%               <XX+YY+ZZ>
%
%   Usage:
%       [x, fval, exitflag, funcCount] = opt_params(p, params, D)  

%     optfun = @(param0) -real(1/2 * (1 - calc_exp(p, param0, D, 2, 2) - ...
%                              calc_exp(p, param0, D, 3, 3) - calc_exp(p, param0, D, 4, 4)));
    optfun = @(param) calc_Heisenberg_exp(p, param, D);
    options = optimoptions('fminunc','GradObj','off','Hessian','off','Display','off',...
        'TolX',1e-5,'TolFun',1e-5, 'Algorithm', 'quasi-newton',...
        'MaxFunEvals', Inf, 'MaxIter', Inf);
    [x, fval, exitflag, output] = fminunc(optfun, params, options);
    funcCount = output.funcCount;
end
