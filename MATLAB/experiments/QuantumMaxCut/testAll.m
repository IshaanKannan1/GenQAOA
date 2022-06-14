% First make graphs and associated random vecs
n = 15;
pmax = 5;
d = 4;
numG = 5;
wGs = {};
As = {};
vs = {};
randns = {};
up = [1;0];
down = [0;1];
minval = intmax;
best = zeros(5,1);
orderlist = [['2','1','3','4'];['4','1','3','2'];['4','1','2','3'];['2','1','4','3']];
seeds = 1:n;

for i = 1:numG
    rng(seeds(i));
    v = normrnd(0, 1, [n,3]);
    v = normr(v);
    [wG, A] = randRegGraph(n,d);
    theta = pi - acos(v(:,3));
    phi = -atan2(v(:,2), v(:,1));  
    psi = cos(theta(1)/2).*down + (cos(phi(1)) + 1i*sin(phi(1)))*sin(theta(1)/2).*up;
    for j = 2:n
        q = cos(theta(j)/2).*down + (cos(phi(j)) + 1i*sin(phi(j)))*sin(theta(j)/2).*up;
        psi = kron(psi, q);
    end
    randns{i} = psi;
    vs{i} = v;
    wGs{i} = wG;
    % As{i} = A;
end

params = nan(pmax, 2, 4, 1);
for p = 1:pmax
    fprintf('p= %i', p);
    for init = 1:2
        if init == 1
            psi0 = ones(2^n, 1)/sqrt(2^n);
        end
        for orders = 1:4
            gtot=0;
            for i = 1:numG
                if init == 2
                    psi0 = randns{i};
                end
                J = [wGs{i}, -ones(size(wGs{i},1))];
                [QAOAhelperfcn, HamObj, HamC, HamZ, HamB, HamD, EvolC, EvolZ, EvolB, EvolD] = SetupQMCHams(n, J, ones(n,1), psi0, vs{i}, [orderlist(orders,:)]);
                for t = 1:6
                    param0 = 2*rand(p,4)-1;
                    myfun = @(param) QAOAhelperfcn(p, param);
                    options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
                    'TolX',1e-5,'TolFun',1e-5, 'Algorithm', 'quasi-newton',...
                    'MaxFunEvals', Inf, 'MaxIter', Inf);
                    [x, fval] = fminunc(myfun, param0, options);

                    if fval < minval
                        minval = fval;
                    end
                end 
                gtot = gtot + minval;
            end
            params(p, init, orders, 1) = gtot/numG;
        end
    end  
end
save('paramsbig.mat', 'params');
