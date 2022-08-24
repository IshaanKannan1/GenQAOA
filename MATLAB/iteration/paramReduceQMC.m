function paramOut = paramReduceQMC(param, D, verbose)
%paramReduceQMC reduces the redundacies in the 4-driver QAOA parameters
%    for quantum MaxCut on a large girth (D+1) regular graph
%
%   Usage: paramOut = paramReduceQMC(param, D)
%          paramOut = paramReduceQMC(param, D, verbose)
%
%   Argument:
%       param = input QAOA parameters in the form of a px4 matrix
%       D = graph branching degree
%       verbose = true if want additional debugging information printed
%                 (optional, default = false)

    if nargin < 3
        verbose = false;
    end


    if size(param,2) ~= 4
        error('Invalid input! param should be px4 matrix')
    end

    delta_tolerance = 1e-4;

    p = size(param,1);

    param = param/pi; % convert it to units of pi

    param = restrictToMiddle(param);

    % restrict alpha to be [-pi/4, pi/4) if D is odd
    if mod(D, 2) == 1
       param(:,1) = mod(param(:,1), 0.5);
       I_fix = param(:,1)>0.25;
       param(I_fix,1) = param(I_fix,1)-0.5;
    else 
        if abs(param(1,1)) > 0.25 % |alpha_1| > pi/4
            if param(1,1) > 0.25
                param(1,1) = param(1,1) - 0.5;
            else
                param(1,1) = param(1,1) + 0.5;
            end
            param(1,2) = -param(1,2);
            param(1,3) = restrictToMiddle(param(1,3) + 0.5);
        end
    end

    % enforce that alpha_1 > 0
    if param(1,1) < 0
        param = -param;
    end

    alphas = param(:,1);
    betas  = param(:,2);
    gammas = param(:,3);
    deltas = param(:,4);

    for ind = 1:p-1

        % try to foce alpha to be in [-pi/4, pi/4]
        if mod(D,2) == 0 && abs(alphas(ind)) > 0.25
            debug_print('fixing alpha_%d\n', ind)
            if alphas(ind) > 0.25
                alphas(ind) = alphas(ind)  - 0.5;
            elseif alphas(ind) < -0.25
                alphas(ind) = alphas(ind)  + 0.5;
            end            
            betas(ind) = -betas(ind);
            gammas(ind) = restrictToMiddle(gammas(ind) - 0.5);
        end

        % check if delta_i = 0 before proceeding
        if deltas(ind) > delta_tolerance
            debug_print('Skipping cuz delta_%d != 0\n', ind)
            continue
        end

        % fix beta_i >= 0
        if betas(ind) < 0
            debug_print('fixing beta_%d\n', ind)
            betas(ind) = betas(ind) + 0.5;
            gammas(ind) = -gammas(ind);
            betas(ind+1) = restrictToMiddle(betas(ind+1) - 0.5);
        end

        % fix gamma_i >= 0
        if gammas(ind) < 0
            debug_print('fixing gamma_%d\n', ind);
            gammas(ind) = gammas(ind) + 0.5;
            if mod(D, 2) == 0
                if alphas(ind+1) > 0.25
                    alphas(ind+1) = alphas(ind+1) - 0.5;
                elseif alphas(ind+1) < -0.25
                    alphas(ind+1) = alphas(ind+1) + 0.5;
                else
                    betas(ind+1) = -betas(ind+1);
                    gammas(ind+1) = restrictToMiddle(gammas(ind+1)-0.5);
                end
            else
                betas(ind+1) = -betas(ind+1);
                gammas(ind+1) = restrictToMiddle(gammas(ind+1)-0.5);
            end
        end
    end

    if mod(D,2) == 0
        if alphas(p) > 0.25
            alphas(p) = alphas(p) - 0.5;
            betas(p) = -betas(p);
            gammas(p) = restrictToMiddle(gammas(p)+0.5);
        elseif alphas(p) < -0.25
            alphas(p) = alphas(p) + 0.5;
            betas(p) = -betas(p);
            gammas(p) = restrictToMiddle(gammas(p)-0.5);
        end
    end

    paramOut = [alphas,betas,gammas,deltas]*pi;
    
    function debug_print(varargin)
        if verbose
            fprintf(varargin{:});
        end
    end
end


function y = restrictToMiddle(x)
%restrict input to be in [-0.5, 0.5] assuming period=1
    y = mod(x,1);
    y(y>0.5) = y(y>0.5) - 1;
end

