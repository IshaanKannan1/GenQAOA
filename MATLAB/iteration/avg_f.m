function fbar = avg_f(z, params, p)
%avg_f Uses a spherical design to average f over all n's, returns a column
%   of fbar(z) values at the given list of bitstrings z and QAOA parameters
%
%   Usage: fbar = avg_f(z, params, p)
%

    path = [fileparts(mfilename('fullpath')) '/SS31-Mar-2016'];
    name = strcat('ss', num2str(2*p+1,'%03.f'));
    fpattern = fullfile(path, strcat(name, '.*'));
    file = dir(fpattern);
    
    pts = load([file(1).folder filesep file(1).name]);
    total = zeros(size(z, 1), 1);
    for j = 1:length(pts)
        total = total + comp_f_iterfun(pts(j, :), z, params, p);
    end
    fbar = total./length(pts);
end
