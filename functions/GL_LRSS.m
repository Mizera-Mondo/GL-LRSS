function [L, X] = GL_LRSS(Y, options)
%GL_LRSS Implementation of GL-LRSS.
% Refer to the paper for detailed explanation for arugments.
arguments
    Y double
    options.R = 0
    options.alpha = 0.1
    options.beta = 0.1
    options.gamma = 0.1
    options.K = 1000
    options.epsilon = 0.5
    options.rho = 0.5
    options.tol = 1e-5
    options.maxIter = 1000
end

%% Preprocess
[a, b] = size(Y);
X = Y;
L = zeros(a);

% Keyboard-friendly names
if options.R == 0
    R = eye(a);
else
    R = options.R;
end

alpha = options.alpha;
beta = options.beta;
gamma = options.gamma;
K = options.K;
epsilon = options.epsilon;
rho = options.rho;
tol = options.tol;
maxIter = options.maxIter;

% Operators
B = zeros(b);
B(2:end, 2:end) = eye(b - 1);
D = @(X) X - R*X*B;


DX = D(X);
CX = DX*DX';


%% Estimation Iteration
iterCount = 1;
isConverge = false;
isMaxIter = false;

while ~isConverge && ~isMaxIter
    L_old = L;
    X_old = X;
    
    tic;
    %% Graph topology refinement
    Z = L;
    Xi = Z - L;
    admmCount = 1;
    isAdmmConverge = false;
    isAdmmMaxIter = false;


    while ~isAdmmConverge && ~isAdmmMaxIter
        L_old_ADMM = L;
        
        % Update of L
        L = (rho*Z + Xi - alpha*CX)/(2*beta + rho);

        % Update of Z
        % FINISHED 0404 TODO: Euclidean projector
        Z = laplacianProjection(L - 1/rho*Xi);

        % Update of Xi
        Xi = Xi + rho*(Z - L);

        % Convergence check
        isAdmmConverge = norm(L_old_ADMM - L,'fro')/norm(L, 'fro') < tol;
        isAdmmMaxIter = admmCount >= maxIter;
        admmCount = admmCount + 1;
    end
    
    toc;

    tic;
    %% Low-rank component estimation
    %'TODO: USING lowRankEstimation
    X = lowRankEstimation(X, D, L, R, B, Y, alpha, gamma);
% %    disp('Low-rank Estimation');
%     P = X;
%     Q = X - P;
%     admmCount = 1;
%     isAdmmConverge = false;
%     isAdmmMaxIter = false;
% 
% 
%     while ~isAdmmConverge && ~isAdmmMaxIter
% %        disp(['Iter: ' num2str(admmCount)]);
%         X_old_ADMM = X;
% 
%         % Update of X
%         count = 1;
% 
%         gradX = @(X) 2*D(X - Y) - 2*R*D(X - Y)*B' + rho*(X - P) ...
%                      + Q + 2*alpha*(L*D(X) - R*L*X*B' + L*X*(B*B'));
%         phi = 2*D(Y) - 2*R*D(Y)*B' +rho*P - Q;
% 
%         % X = zeros(size(X));
%         gX = gradX(X);
%         dX = -gX;
% 
%         % Optimize X with Flecther-Reeves conjugate gradient method
%         while count < 10
% 
%             mu = -1*trace(dX'*gX)/trace(dX'*(gX + phi)); % The minus sign of mu may be a mistake
%             X = X + mu*dX;
% 
%             gXnew = gradX(X);
%             theta = norm(gXnew, 'fro')/norm(gX, 'fro');
%             theta = theta^2;
% 
%             gX = gXnew;
%             dX = -gX + theta*dX;
%             count = count + 1;
% 
%         end
% 
% 
%         % Update of auxiliary P
%         P = singularValueThreshold(X + 1/rho*Q, gamma/rho);
% 
%         % Update of Q
%         Q = Q + rho*(X - P);
% 
% 
%         % Convergence check
%         isAdmmConverge = norm(X_old_ADMM - X,'fro')/norm(X, 'fro') < tol;
%         isAdmmMaxIter = admmCount >= maxIter;
%         admmCount = admmCount + 1;
%     end    
    
    toc;
    %% Terminate condition check
    % isConverge = (norm(L_old - L, 'fro')/norm(L, 'fro') < tol) && ...
    %             (norm(X_old - X, 'fro')/norm(X, 'fro') < tol);
    isConverge = (norm(L_old - L, 'fro')/norm(L, 'fro') < tol);

    isMaxIter = iterCount >= maxIter;
    iterCount = iterCount + 1;
    
    disp(num2str(norm(L_old - L, 'fro')/norm(L, 'fro')));
end
end
