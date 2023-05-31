% Noiseless signal constrcution
n = 20;
m = 1000;
R = eye(n) - diag(0.05*rand(n, 1));
A = rand_ugraph(n, ceil(n^2/10), 0.1, 0.1);
L = diag(A*ones(n, 1)) - A;
L = L/trace(L)*n;

[vec, val] = eig(L);
[vec, val] = sortEigen(vec, val, 'ascend');
U = vec(:, 1:ceil(0.25*n));
sigma = val(1:ceil(0.25*n), 1:ceil(0.25*n));
Z = randn(m, ceil(0.25*n))*chol(sigma+0.001*eye(ceil(0.25*n)));
Z = Z';
V = U*Z;

X(:, 1) = V(:, 1);
for i = 2:m
    X(:, i) = R*X(:, i - 1) + V(:, i);
end

Lest = GLLRSS_GRONLY(X, R = R, beta = 0.05, tol = 1e-7);
close all;
figure; imagesc(L); colorbar; title('Ground Truth');
figure; imagesc(Lest); colorbar; title('Estimated');

% save('temp.mat');