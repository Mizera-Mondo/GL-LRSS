function currentValue = lineSearchArminjo(startValue, currentGrad, targetFun, adoptCoff, maxIter)
%LINESEARCH Line searcher with arminjo criteria.
i = 1;
isMaxIterReach = false;
isLSConverge = false;
while ~isMaxIterReach && ~isLSConverge
    deltaV = -1*(0.9)^i*currentGrad;
    currentValue = startValue + deltaV;
    isLSConverge = ((targetFun(currentValue) - targetFun(startValue))) <= (adoptCoff*(trace(deltaV'*currentGrad)));
    isMaxIterReach = i >= maxIter;
    i = i + 1;
end
end
