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
nodeNum = 20;
usedEigNum = 15;
signalLength = 2000;
noiseCov = 0.01;
rPertubation = 0.01;

[Y, A, R] = genRandomSignal(nodeNum, usedEigNum, signalLength, noiseCov, rPertubation);

L = diag(sum(A)) - A;
B = zeros(signalLength);
B(2:end, 2:end) = eye(signalLength - 1);
D = @(X) X - R*X*B;
gamma = 0.001;
countAlpha = 1;
countBeta = 1;
errorGridAlpha = [];
errorGridBeta = [];
errorRec = [];
for alpha = 0.01:0.3:1
    for beta = 0.1:0.5:5
        [Lest, X] = GL_LRSS(Y, R = R, beta = beta, gamma = gamma, tol = 1e-4);
        errorGridAlpha(countAlpha, countBeta) = alpha;
        errorGridBeta(countAlpha, countBeta) = beta;
        errorRec(countAlpha, countBeta) = norm(Lest - L, 'fro');
        countBeta = countBeta + 1;
    end
    countBeta = 1;
    countAlpha = countAlpha + 1;
end