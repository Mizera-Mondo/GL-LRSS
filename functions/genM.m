function M = genM(DX)
%GENM Mij = 1/2||d_i - d_j||_2^2 
    [nodeNum, ~] = size(DX);
    M1 = DX*DX';
    M2 = repmat(diag(M1), 1, nodeNum);
    M = (M2 + M2' - 2*M1)./2;
end

