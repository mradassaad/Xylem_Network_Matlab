function pg = vizCav(cavGraph)
LWidths = 3*cavGraph.Edges.Weight/max(cavGraph.Edges.Weight);
%Color inlet conduits in green and outlet conduits in red
hIn = find(cavGraph.Nodes.isInlet);
hOut = find(cavGraph.Nodes.isOutlet);
hRedundant = find(cavGraph.Nodes.isRedundant);
hEmbolized = find(cavGraph.Nodes.isEmbolized);
figure;
pg = plot(cavGraph, 'layout', 'layered',...
    'Sources', hIn, 'Sinks',  hOut, 'Direction', 'up',...
    'LineWidth', LWidths, 'MarkerSize', 5);
highlight(pg, hIn, 'NodeColor', 'g');
highlight(pg, hOut, 'NodeColor', 'r');
highlight(pg, hRedundant, 'NodeColor', 'k');
highlight(pg, hEmbolized, 'NodeColor', 'm');
end