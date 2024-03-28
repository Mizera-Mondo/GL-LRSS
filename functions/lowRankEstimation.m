function [X, P] = lowRankEstimation(X, D, L, R, B, Y, alpha, gamma, options)
%LOWRANKESTIMATION Solve ||D(X - Y)||_F^2 + alpha(D(X)'*L*D(X)) +
%gamma||P||_*, X = P
%   Using ADMM
arguments
    X, D, L, R, B, Y, alpha, gamma double
    options.tol = 1e-4
    options.method = "GD"
    options.maxIter = 1000
    options.rho = 1

end
rho = options.rho;
P = X;
Q = X - P;
gradX = @(X, P, Q) 2*D(X - Y) - 2*R*(D(X - Y))*B' + rho*(X - P) + Q + ...
    2*alpha*(L*D(X) - R*L*X*B' + L*X*(B*B'));
targetFunX = @(X, P, Q) (norm(D(X - Y), 'fro'))^2 + alpha*trace((D(X))'*L*D(X)) + ...
    rho/2*(norm(X - P + Q/rho, 'fro'))^2;

iter = 1;
isADMMConverge = false;
isMaxIter = false;

%% Estimation Iterations
while ~isADMMConverge && ~isMaxIter
    X_old = X;
    P_old = P;
% Update of X
    targetFunX_ = @(X) targetFunX(X, P, Q);
    isXConverge = false;
    isXMaxIter = false;
    xIter = 1;
    while ~isXConverge && ~isXMaxIter
        X_sub_old = X;
        X = lineSearchArminjo(X, gradX(X, P, Q), targetFunX_, 1e-3, 100);
        isXConverge = norm(X - X_sub_old, 'fro')/norm(X_sub_old, 'fro') <= options.tol;
        isXMaxIter = xIter >= options.maxIter;
        xIter = xIter + 1;
    end
% Update of P
    P = singularValueThreshold(X + Q/rho, gamma/rho);
% Update of Q
    Q = Q + rho*(X - P);
% Terminate Condition Check
    changeRateX = norm(X - X_old, 'fro')/norm(X_old, 'fro');
    changeRateP = norm(P - P_old, 'fro')/norm(P_old, 'fro');
    isADMMConverge = changeRateX <= options.tol && changeRateP <= options.tol;
    isMaxIter = iter >= options.maxIter;
    iter = iter + 1;
end

end