% Noiseless signal constrcution
% n = 20;
% m = 1000;
% R = eye(n) - diag(0.05*rand(n, 1));
% A = rand_ugraph(n, ceil(n^2/10), 0.1, 0.1);
% L = diag(A*ones(n, 1)) - A;
% L = L/trace(L)*n;
% 
% [vec, val] = eig(L);
% [vec, val] = sortEigen(vec, val, 'ascend');
% U = vec(:, 1:ceil(0.25*n));
% sigma = val(1:ceil(0.25*n), 1:ceil(0.25*n));
% Z = randn(m, ceil(0.25*n))*chol(sigma+0.001*eye(ceil(0.25*n)));
% Z = Z';
% V = U*Z;
% 
% X(:, 1) = V(:, 1);
% for i = 2:m
%     X(:, i) = R*X(:, i - 1) + V(:, i);
% end
nodeNum = 30;
usedEigNum = 5;
signalLength = 100;
noiseCov = 1;
rPertubation = 0.01;

[Y, A, R] = genRandomSignal(nodeNum, usedEigNum, signalLength, noiseCov, rPertubation);

B = zeros(signalLength);
B(2:end, 2:end) = eye(signalLength - 1);
D = @(X) X - R*X*B;
alpha = 0.1;
beta = 0.1;
gamma = 0.01;
targetFunction = @(L, X) (norm(D(X - Y), 'fro'))^2 + alpha*trace((D(X))'*L*D(X)) + beta*(norm(L, 'fro'))^2 + gamma*(nuclearNorm(X));


[Lest, X] = GL_LRSS(Y, R = R, beta = beta, tol = 1e-4);
close all;
figure; imagesc(L); colorbar; title('Ground Truth');
figure; imagesc(Lest); colorbar; title('Estimated');

% save('temp.mat');