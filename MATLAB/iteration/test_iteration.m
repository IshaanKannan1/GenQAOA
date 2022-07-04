% set up graph and random vectors
p = 2;
D = 2;
k = 2;
n = 2*(2^(k+1) - 1);
[wg, A] = gen_bbst(k + 1);
v = normrnd(0, 1, [n,3]);
v = normr(v);

% set up initial state
up = [1;0];
down = [0;1];
theta = pi - acos(v(:,3));
phi = -atan2(v(:,2), v(:,1));  
psi = cos(theta(1)/2).*down + (cos(phi(1)) + 1i*sin(phi(1)))*sin(theta(1)/2).*up;
for j = 2:n
    q = cos(theta(j)/2).*down + (cos(phi(j)) + 1i*sin(phi(j)))*sin(theta(j)/2).*up;
    psi = kron(psi, q);
end

J = [wg, -ones(size(wg,1), 1)];

[QAOAhelperfcn, HamObj, HamC, HamZ, HamB, HamD, EvolC, EvolZ, EvolB, EvolD] = SetupQMCHams(n, J, ones(n,1), psi, v, 1);
param0 = ones(p,4);
% myfun = @(param) QAOAhelperfcn(p, param);
% options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
%     'TolX',1e-5,'TolFun',1e-5, 'Algorithm', 'quasi-newton','PlotFcns',{@optimplotfval, @optimplotstepsize},...
%     'MaxFunEvals', Inf, 'MaxIter', Inf);
% [x, fval] = fminunc(myfun, param0, options);
[F, Fgrad, psi] = QAOAhelperfcn(p, param0);
iter_val = 1/2 * (1 - calc_exp(p, param0, D, 2, 2) - ...
            calc_exp(p, param0, D, 3, 3) - calc_exp(p, param0, D, 4, 4));
-F
iter_val
% z = [1 1 1 1 1 1 1 1];
% zm = -z;
% avg_f(z, x(:,2:end), p)
% avg_f(zm, x(:,2:end), p)
