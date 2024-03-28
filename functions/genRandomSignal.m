function [x, R] = genRandomSignal(L, usedEigNum, signalLength, noiseCov, rPerturbation)
    
    % generate R
    [nodeNum, ~] = size(L);
    R = eye(nodeNum) + diag(-2*rPerturbation*ones(nodeNum, 1));
    
    % generate vt
    [vec, val] = eig(L);
    [vec, val] = sortEigen(vec, val, 'ascend');
    vecT = vec(:, 1:usedEigNum);
    valT = val(1:usedEigNum, 1:usedEigNum);
    covZ = pinv(neatZero(valT));


    z = randn(usedEigNum, signalLength);
    z = sqrt(covZ)*z;
    v = vecT*z;
    
    % signal construction
    x = zeros(nodeNum, signalLength);
    x(:, 1) = v(:, 1);
    for t = 2:signalLength
        x(:, t) = R*x(:, t - 1) + v(:, t);
    end

    % adding white noise
    x = x + sqrt(noiseCov)*eye(nodeNum)*randn(nodeNum, signalLength);
end