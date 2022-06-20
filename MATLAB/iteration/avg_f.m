function avg = avg_f(z, params, p)

    path = '/Users/ishaan/Documents/MATLAB/GenQAOA-master/MATLAB/iteration/SS31-Mar-2016';
    name = strcat('ss', num2str(2*p+1,'%03.f'));
    fpattern = fullfile(path, strcat(name, '.*');
    file = dir(fpattern);
    pts = load(file(1).name);
    total = 0;
    for j = 1:length(pts) %fine since p > 0 so length >= 3
        total = total + comp_f_iterfun(pts(j, :), z, params, p);
    end
    avg = total/length(pts);

end