function paramspace(mydata, p, D)
    
    energy = reshape([mydata.energy], D, []);
    Ds = reshape([mydata.D], D, []);
    Ds = Ds(:,1);
    for i = 1:length(Ds)
        subplot(3, 3, i);
        bestpars = mod(get_best_pars(mydata, D, i), 2*pi)-pi;
        a1 = bestpars(1, 1:4:end);
        b1 =  bestpars(1, 2:4:end);
        g1 = bestpars(1, 3:4:end);
        d1 = bestpars(1, 4:4:end);
        plot(1:4, [a1;b1;g1;d1], '-o');
        title(sprintf('D = %d', Ds(i)));
        set(gca, 'xscale', 'log');
    end
end