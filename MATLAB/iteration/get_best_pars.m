function [pars, small] = get_best_pars(mydata, Dtot, ind)
    energy = reshape([mydata.energy], Dtot, []);
    temp = sort(energy(ind, :));
    pars = [mydata(ind, [mydata(ind,:).energyda] <= temp(1) + 1e-5).finalpars];
    small = [];
    for i = 1:size(pars, 2)/4
        if pars(:, 4*(i-1) + 1) < pi/4 & pars(:, 4*(i-1) + 1) > 0
            small = [small, pars(:, 4*(i-1) + 1:4*(i-1)+4)];
        end
    end
end