function pg = vizCond(condGraph)
ggCond = condGraph;
LWidths = 3*ggCond.Edges.Weight/max(ggCond.Edges.Weight)+1e-20;
hRedundant = find(condGraph.Nodes.isRedundant);
hEmbolized = find(condGraph.Nodes.isEmbolized);
figure;
pg = plot(ggCond, 'XData', ggCond.Nodes.XData,...
    'YData', ggCond.Nodes.YData, 'ZData', ggCond.Nodes.ZData,...
    'MarkerSize',.1,'LineWidth',LWidths, 'MarkerSize', 2);
highlight(pg, hRedundant, 'NodeColor', 'k')
highlight(pg, hEmbolized, 'NodeColor', 'm');
end