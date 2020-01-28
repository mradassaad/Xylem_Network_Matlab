function coeff = clusCoeff(g)
    nodes = g.Nodes;
    connPair = zeros(1,height(nodes));
    for i = 1 : height(nodes)
        neighbs = neighbors(g,i);
        while length(neighbs) > 1
            n =  neighbs(1);
            %compute number of connected pairs
            connPair(i) = connPair(i) +...
                sum(ismember(neighbors(g,n),neighbs(2:end)));
            neighbs = neighbs(2:end);
        end
    end
    degs = degree(g);
    coeff = 2*sum(connPair)/sum(degs.*(degs-1));
end