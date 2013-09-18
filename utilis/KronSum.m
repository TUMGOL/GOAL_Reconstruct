function c = KronSum( a, b )
c = kron( a, ones(size(b)) ) + kron( ones(size(a)), b );