function data = find_opt_pars(p, nTrials, Ds, runID, varargin)
% Uses gradient descent on num_p randomly sampled points from a 4p-dimensional cube of side length
% pi from -pi/2 to pi/2 to compute optimal parameters for QMC QAOA iteration for given degrees.
% Saves resulting data in a struct file.
% Arguments:
%   - p: QAOA depth
%   - nTrials: Desired number of initial points
%   - D: List of positive integers indicating degrees for which we want to optimize
%
% Returns:
%   best_v, best_pars, best_D

    fname = [fileparts(mfilename('fullpath')), ...
             sprintf('/paramsearch/p=%d_ID=%d.mat',p,runID)];
%     Ds = [1, 2, 3, 5, 10, 20, 100];

    rng(20220723+997*runID+p);

    fprintf('Saving to %s\n', fname);

    data(length(Ds), nTrials).D = [];
    data(length(Ds), nTrials).energy = [];
    data(length(Ds), nTrials).exitflag = [];
    data(length(Ds), nTrials).initpars = [];
    data(length(Ds), nTrials).finalpars = [];
    data(length(Ds), nTrials).funcCount = [];
    data(length(Ds), nTrials).realTime = [];

    if nargin > 4
        varargin = varargin{1};
        nTrials = size(cell2mat(varargin(1)), 2) / 4;
    end

    for i = 1:length(Ds)
        for j = 1:nTrials
            tic;
%             disp([i,j])
            if nargin > 4
                param = cell2mat(varargin(i));
                param0 = param(:, 4*j-3: 4*j);
                param0(param0 == 0) = (2*rand(1) - 1) * .02
            else
                param0 = (2*rand(p, 4)-1) .* pi/2;
            end
            [x, fval, exitflag, funcCount] = opt_params(p, param0, Ds(i));
            data(i, j).D = Ds(i);
            data(i, j).energy = fval;
            data(i, j).exitflag = exitflag;
            data(i, j).initpars = param0;
            data(i, j).finalpars = x;
            data(i, j).funcCount = funcCount;
            data(i, j).realTime = toc;
            fprintf('p=%d, D=%d, trial=%d, took %0.2f s\n', p, Ds(i), j, toc);
        end
    end
    save(fname, 'data');
%     [best_v, I] = min([data.energy], [], 2);
%     best_v = -best_v;
%     best_pars = data(I).finalpars;
%     best_D = data(I).D;
end
