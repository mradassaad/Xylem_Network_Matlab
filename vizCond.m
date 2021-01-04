function pg = vizCond(condGraph)
ggCond = condGraph;
LWidths = 3*ggCond.Edges.Weight/max(ggCond.Edges.Weight)+1e-20;
% hRedundant = find(condGraph.Nodes.isRedundant);
% hEmbolized = find(condGraph.Nodes.isEmbolized);
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
% highlight(pg, hRedundant, 'NodeColor', 'k')
% highlight(pg, hEmbolized, 'NodeColor', 'm');
iccBool = ggCond.Edges.Type == "ICC";
s = ggCond.Edges.EndNodes(iccBool,1);
t = ggCond.Edges.EndNodes(iccBool,2);
highlight(pg, s, t, 'EdgeColor', 'red')

ceBool = ggCond.Edges.Type == "CE";
s = ggCond.Edges.EndNodes(ceBool,1);
t = ggCond.Edges.EndNodes(ceBool,2);
highlight(pg, s, t, 'EdgeColor', 'black')
end