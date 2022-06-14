% use QAOA script function to test on an ensemble of different drivers and
% try different orders. Should edit test_func to take in different driver
% order. 

N = 10;
numGs = 50;
deg = 3;
[wG, A] = randRegGraph(N, deg);
wGs = [wG];
As = [A];
for i = 1:numGs
    [a,b] = randRegGraph(N,deg);
    wGs = [wGs, a];
    As = [As, b];
end

avg = 0;
v = perms(['1','2','3','4']);
for j = 1:length(v)
    avg = 0;
    for i = 1:numGs
        avg = avg + QMC_test_func(wGs(:, 2*i-1:2*i), As(:,1 + (i-1)*N:N*i), v(j,:));
    end
    avg = avg / numGs;
    fprintf('Average over %d %d-regular random graphs: %d, Order: %s',numGs, deg, avg, v(j,:));
    disp('*****');
end


    