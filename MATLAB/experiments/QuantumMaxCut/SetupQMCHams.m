function [QAOAhelperfcn, HamObj, ...
    HamC, HamZ, HamB, HamD, EvolC, EvolZ, EvolB, EvolD] =  SetupQMCHams(N, J, h, psi0, v, oneEdge, varargin)
%SetupQMCHams sets up the Hamiltonians and evolution functions for QAOA
%           in one-dimension using ZZ, Z, and X-type driver to optimize the
%           Quantum MaxCut Hamiltonian problem
%
%   Arguments: 
%       N = system size
%       J = an |E|x3 array of edges and weights whos rows are [i, j, J_ij]
%       h = a |V|x1 array specifying coefficient in the Z driver
%                   so that HamZ = \sum_i h_i Z_i as driver,
%            (default or [] (empty) if using \sum_i Z_i as driver)
%
%   Returns:
%       QAOAhelperfcn = @(p, param) an evolution function based on
%                       MultiQAOAGrad, returns <HamObj> and its gradient
%                       for a given p and a set of parameters param
%       HamObj = the XXZ objective Hamiltonian, which is
%                \sum_ij J_ij/2 * (1 - XX - YY - ZZ)_{i,j}
%       HamC = ZZ-type driver (\sum_ij J_ij Z_i Z_j)
%       HamZ = Z-type driver (\sum_i h_i Z_i)
%       HamB = X-type driver (\sum_i X_i)
%       EvolC = @(psi, gamma) evolves a given psi by exp(-1i*gamma*HamC)
%       EvolZ = @(psi, alpha) evolves a given psi by exp(-1i*alpha*HamZ)
%       EvolB = @(psi, beta) evolves a given psi by exp(-1i*beta*HamB)
%


%% Set up Hamiltonian
sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));

% objective Hamiltonian
HamObj = sparse(0);

% For testing iteration, only care about edge between middle vertices,
% which will always be in the same place, (N/2)
if oneEdge == 1
    Ia = J(N/2, 1);
    Ib = J(N/2, 2);
    HamObj = HamObj + J(ceil(N/2),3)/2 * (speye(2^N) - Ham2LTerm(sx, sx, Ia, Ib, N) -...
            Ham2LTerm(sy, sy, Ia, Ib, N) - Ham2LTerm(sz, sz, Ia, Ib, N));
else
    for ind = 1:size(J, 1)
        Ia = J(ind,1);
        Ib = J(ind, 2);
        HamObj = HamObj + J(ind,3)/2 * (speye(2^N) - Ham2LTerm(sx, sx, Ia, Ib, N) -...
            Ham2LTerm(sy, sy, Ia, Ib, N) - Ham2LTerm(sz, sz, Ia, Ib, N));
    end
end


% ZZ driver

HamC = sparse(0);

for ind = 1:size(J,1)
    Ia = J(ind,1);
    Ib = J(ind,2);
    HamC = HamC + J(ind, 3) * Ham2LTerm(sz, sz, Ia, Ib, N);
end

if nargin <= 2 || isempty(h)
    HamZ = krondist(sz, N);
else
    HamZ = krondist(sz, N, h);
end

% random pauli driver
% v = normrnd(0, 1, [N,3]);
% v = v/norm(v);
HamD = krondist(sx, N, v(:,1)) + krondist(sy, N, v(:,2)) + krondist(sz, N, v(:,3));
%

HamB = krondist(sx, N);
HamZdiag = full(diag(HamZ));
HamCdiag = full(diag(HamC));

EvolB = @(psi, beta) EvolHamB(N, beta, psi, false);
EvolC = @(psi, gamma) exp(-1i*gamma*HamCdiag).*psi;
EvolZ = @(psi, alpha) exp(-1i*alpha*HamZdiag).*psi;
EvolD = @(psi, delta) PauliRotations(N, delta, v, psi);


% psi0 = ones(2^N,1)/sqrt(2^N);

if nargin > 6
    %varargin order: X, ZZ, Z, RP (BACD)
    x = num2cell(cell2mat(varargin));
    M = containers.Map(x, {HamB, HamC, HamZ, HamD});
    E = containers.Map(x, {EvolB, EvolC, EvolZ, EvolD});    
    Hams = {M('1'), M('2'), M('3'), M('4')};
    Evols = {E('1'), E('2'), E('3'), E('4')};
    QAOAhelperfcn = @(p, param) MultiQAOAGrad(p, HamObj, Hams, param, psi0, Evols);
else 
    QAOAhelperfcn = @(p, param) MultiQAOAGrad(p, HamObj, {HamC, HamB, HamZ, HamD}, param, psi0, {EvolC, EvolB, EvolZ, EvolD});
end
