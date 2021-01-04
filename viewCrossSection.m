function viewCrossSection(xn,row)
%VIEWCROSSSECTION Allows to see 2D cross-section
%   
ggCond = rmnode(xn.gCond, find(xn.gCond.Nodes.XData ~= row));

LWidths = 3*ggCond.Edges.Weight/max(ggCond.Edges.Weight)+1e-20;
% hRedundant = find(condGraph.Nodes.isRedundant);
% hEmbolized = find(condGraph.Nodes.isEmbolized);
figure;
r0 = 10;
r = r0 + ggCond.Nodes.YData;
thetaStep = (2*pi)/max(ggCond.Nodes.ZData);
thetaVal = (0:thetaStep:2*pi - thetaStep);
theta = thetaVal(ggCond.Nodes.ZData)';

xNew = r.*cos(theta);
yNew = r.*sin(theta);
pg = plot(ggCond, 'XData', xNew,...
    'YData', yNew,...
    'MarkerSize',2,'LineWidth',LWidths,'NodeLabel',{});

iccBool = ggCond.Edges.Type == "ICC";
s = ggCond.Edges.EndNodes(iccBool,1);
t = ggCond.Edges.EndNodes(iccBool,2);
highlight(pg, s, t, 'EdgeColor', 'red')

maxComp = max(ggCond.Nodes.component);
minComp = min(ggCond.Nodes.component);
colors = {[0, 0.4470, 0.7410], [0, 0, 1], [0.8500, 0.3250, 0.0980],...
    [0, 0.5, 0], [0.9290, 0.6940, 0.1250], [1, 0, 0], [0.4940, 0.1840, 0.5560],...
    [0, 0.75, 0.75], [0.4660, 0.6740, 0.1880], [0.75, 0, 0.75],...
    [0.3010, 0.7450, 0.9330], [0.75, 0.75, 0], [0.6350, 0.0780, 0.1840],[0.25, 0.25, 0.25]};
colors = repmat(colors,1, 6);
for i = minComp : maxComp
    
    highlight(pg, ggCond.Nodes.component==i, 'NodeColor',...
        colors{i - minComp + 1});
end

end

