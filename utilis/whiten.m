function [X] = whiten(X,fudgefactor)
X = bsxfun(@minus, X, mean(X));
A = X'*X;
[V,D] = eig(A);
X = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
if numel(X(isnan(X)))
    hallo=1;
end
X = X(:);
end