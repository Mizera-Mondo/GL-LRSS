function [x, A, R] = genRandomSignal(nodeNum, usedEigNum, signalLength, noiseCov, rPerturbation)
    A = rand_ugraph(nodeNum, floor(nodeNum*log(nodeNum)/1.5), 0.1, 0.1); % log(nodeNum)/1.5 links per node
    A = A./sum(A, "all")*nodeNum;
    L = diag(sum(A)) - A;
    R = eye(nodeNum) + diag(-2*rPerturbation*ones(nodeNum, 1));
    
    [vec, val] = eig(L);
    [vec, val] = sortEigen(vec, val, 'ascend');
    vecT = vec(:, 1:usedEigNum);
    valT = val(1:usedEigNum, 1:usedEigNum);
    covZ = pinv(neatZero(valT));
    % 'Is the cov right? no inversion?
    % It is right. Pay attention to the zero in valT.

    z = randn(usedEigNum, signalLength);
    z = sqrt(covZ)*z;
    v = vecT*z;

    x = zeros(nodeNum, signalLength);
    x(:, 1) = v(:, 1);
    for t = 2:signalLength
        x(:, t) = R*x(:, t - 1) + v(:, t);
    end
    x = x + sqrt(noiseCov)*eye(nodeNum)*randn(nodeNum, signalLength);
end