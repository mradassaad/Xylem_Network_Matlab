function pg = vizCond(condGraph)
ggCond = condGraph;
LWidths = 3*ggCond.Edges.Weight/max(ggCond.Edges.Weight)+1e-20;
hRedundant = find(condGraph.Nodes.isRedundant);
hEmbolized = find(condGraph.Nodes.isEmbolized);
figure;
zCyl = ggCond.Nodes.XData;
r0 = 10;
r = r0 + ggCond.Nodes.YData;
thetaStep = (2*pi)/max(ggCond.Nodes.ZData);
thetaVal = (0:thetaStep:2*pi - thetaStep);
theta = thetaVal(ggCond.Nodes.ZData)';

xNew = r.*cos(theta);
yNew = r.*sin(theta);
zNew = zCyl;
pg = plot(ggCond, 'XData', xNew,...
    'YData', yNew, 'ZData', zNew,...
    'MarkerSize',.1,'LineWidth',LWidths, 'MarkerSize', 0.1);
highlight(pg, hRedundant, 'NodeColor', 'k')
highlight(pg, hEmbolized, 'NodeColor', 'm');
end