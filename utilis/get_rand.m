function series = get_rand(range,nbr)
    randn('seed',0);
    rand('seed',0);
    series = randperm(range);
    if nargin == 2
        series = series(1:nbr)';
    end
end