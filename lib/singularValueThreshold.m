function M = singularValueThreshold(M, tau)
%SINGULARVALUETHRESHOLD 
[U, Z, V] = svd(M);
Z = sign(Z).*...
    ((abs(Z) - tau) > 0).* ...
    (abs(Z) - tau);
M = U*Z*V';
end
