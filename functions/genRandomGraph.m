function [A, L] = genRandomGraph(nodeNum, options)
    arguments
        nodeNum double
        options.nodeNum = floor(nodeNum*log(nodeNum)/1.5); % log(nodeNum)/2 links per node
    end
    A = rand_ugraph(nodeNum, options.nodeNum, 0.1, 0.1); 
    A = A./sum(A, "all")*nodeNum;
    L = diag(sum(A)) - A;
end