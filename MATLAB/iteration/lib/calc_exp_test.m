p = 6;
params = rand(p,4);

D = 2;
sL = 4;
sR = 4;

%% new version

tic;
A = params(:,1);
A = [A;0;0;-flip(A)];

bList = logical(de2bi(0:2^(2*p)-1, 2*p));
bList = [bList(:,1:p), false(2^(2*p),2), bList(:,p+1:end)];
bList2 = bList; bList2(:,p+1) = true;
bList3 = bList; bList3(:,p+2) = true;
bList4 = bList2; bList4(:,p+2) = true;
zList = 2*[bList;bList2;bList3;bList4] - 1;

H_p = iterate(2*bList-1, params, p, D);
H_p = [H_p;H_p; H_p; H_p];

clear bList bList2 bList3 bList4

fbar = avg_f(zList, params(:, 2:end), p);
hl = comp_h_iterfun(zList, sL, p);
hr = comp_h_iterfun(zList, sR, p);

fprintf('Prep time = %0.6f s\n', toc);

%%
tic;
exp_val = 0;
batch_size = 2^p;
for ind = 1:batch_size:length(zList)
    iL = ind-1 + (1:batch_size);

    exp_val = exp_val + sum( H_p(iL) .* fbar(iL) .* hl(iL) ...
         .* (exp(-1i * (zList(iL, :).* A.') * zList.' ) * (H_p .* fbar .* hr)) );
end

fprintf('Final value time = %0.6f s\n', toc);

%% old version

tic;

A = params(:,1);
A = [A;0;0;-flip(A)];

zList = 2* de2bi(0:2^(2*p+2)-1) -1 ;
H_p = iterate(zList, params, p, D);

fbar = avg_f(zList, params(:, 2:end), p);
hl = comp_h_iterfun(zList, sL, p);
hr = comp_h_iterfun(zList, sR, p);

fprintf('OLD Version: Prep time = %0.6f s\n', toc);

%%
tic;
    old_val = 0;
    for iL = 1:length(zList)
        old_val = old_val + H_p(iL) * fbar(iL) * hl(iL) ...
             * (exp(-1i * (zList(iL, :).*zList) * A).' * (H_p .* fbar .* hr));
    end
fprintf('OLD Version: Final value time = %0.6f s\n', toc);

