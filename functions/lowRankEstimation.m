function [X, P] = lowRankEstimation(X, D, L, R, B, Y, alpha, gamma, options)
%LOWRANKESTIMATION Solve ||D(X - Y)||_F^2 + alpha(D(X)'*L*D(X)) +
%gamma||P||_*, X = P
%   Using ADMM
arguments
    X, D, L, R, B, Y, alpha, gamma double
    options.method = "GD"

end
P = X;
Q = X - P;
rho = 0.5;
gradX = @(X, P, Q) 2*D(X - Y) - 2*R*(D(X - Y))*B' + rho*(X - P) + Q + ...
    2*alpha*(L*D(X) - R*L*X*B' + L*X*(B*B'));
targetFunX = @(X, P, Q) (norm(D(X - Y), 'fro'))^2 + alpha*trace((D(X))'*L*D(X)) + ...
    rho/2*(norm(X - P + Q/rho, 'fro'))^2;

tol = 1e-6;
iter = 1;
maxIter = 1000;
isADMMConverge = false;
isMaxIter = false;

%% Estimation Iterations
while ~isADMMConverge && ~isMaxIter
    X_old = X;
    P_old = P;
% Update of X
    if strcmp(options.method, "GD")
        targetFunX_ = @(X) targetFunX(X, P, Q);
        X = lineSearchArminjo(X, gradX(X, P, Q), targetFunX_, 0.05, 1000);
    end
% Update of P
    P = singularValueThreshold(X + Q/rho, gamma/rho);
% Update of Q
    Q = Q + rho*(X - P);
% Terminate Condition Check
    isADMMConverge = (norm(X - X_old, 'fro'))^2/(norm(X_old, 'fro'))^2 < tol &&...
        (norm(P - P_old, 'fro'))^2/(norm(P_old, 'fro'))^2 < tol;
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end

end