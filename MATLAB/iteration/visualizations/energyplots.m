function energyplots(p, mydata)
    Ds = reshape([mydata.D], 4, []);
    Ds = Ds(:,1);
    energy = reshape([mydata.energy], 4, []);
    
    figure(11);
    for ind = 1:length(Ds)
        subplot(2,2,ind);
        temp = sort(energy(ind,:));
        nCounts = nnz(~isnan(temp));
        plot(temp,'o');
        set(gca,'xscale','log');
        grid on
        xlabel('index of LM (sorted)')
        ylabel('LM value')
        title(sprintf('D=%d,  #min=%d/%d', Ds(ind), nnz(temp <= temp(1) + 1e-5), nCounts));
    end
    
    figure(12);
    % clf
    bestEnergy = min(energy,[],2);
    plot(Ds, -bestEnergy, 'o:b');
    hold on
    plot(Ds, 1./sqrt(Ds), 'r--')
    hold off, grid on
    set(gca,'xscale','log','yscale','log')
    title(sprintf('p=%d', p));
    xlabel('branching degree D')
    ylabel('-\langleXX + YY + ZZ\rangle')
end