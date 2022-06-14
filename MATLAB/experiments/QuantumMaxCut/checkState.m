n = 4;
d = 2;
v = normrnd(0, 1, [n,3]);
v = normr(v);
% for i = 1:n
%     v(i, :) = v(i,:)/norm(v(i,:));
% end
[wG, A] = randRegGraph(n,d);
theta = pi - acos(v(:,3));
phi = -atan2(v(:,2), v(:,1));
up = [1;0];
down = [0;1];

psi = cos(theta(1)/2).*down + (cos(phi(1)) + 1i*sin(phi(1)))*sin(theta(1)/2).*up;
for j = 2:n
    q = cos(theta(j)/2).*down + (cos(phi(j)) + 1i*sin(phi(j)))*sin(theta(j)/2).*up;
    psi = kron(psi, q);
end

sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));
HamD = krondist(sx, n, v(:,1)) + krondist(sy, n, v(:,2)) + krondist(sz, n, v(:,3));

[V, D] = eigs(HamD, 2);
size(HamD)
V(:,2)
-1*psi
a= V(:,2) + psi

norm(V(:,2)-psi)/norm(V(:,2))
norm(V(:,1)-psi)/norm(V(:,1))
