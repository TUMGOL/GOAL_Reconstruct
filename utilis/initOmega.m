function Omega = initOmega(X,rows)
Omega = zeros(rows,size(X,1));
d = size(X,1)-1;
for i=1:rows
    sel = randperm(length(X));
    Omega(i,:) = null(X(:,sel(1:d))')';
end