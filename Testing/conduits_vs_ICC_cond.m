sim = Agg;
condK = zeros(1, length(sim.Conduits));
ICCK = zeros(1, length(sim.Conduits));
for i = 1 : length(sim.Conduits)
    condK(i) = sim.Conduits(i).Kce * sim.Conduits(i).Lce;
    ICCK(i) = sum([sim.Conduits(i).ICConnections.Km]) *...
        sum([sim.Conduits(i).ICConnections.Am]) / sim.Conduits(i).Length;
end

plot(condK, ICCK, '.')
hold on
plot([0 max([condK ICCK])],[0 max([condK ICCK])])

mu = 1.002e-3; %Pa.s
Ds =  {sim.Dcross(1,:), sim.Dcross(2,:),...
    sim.Dcross(3,:), sim.Dcross(4,:)};
As = cell(1, 4);
K_HP = cell(1,4);
K_HP_sum = zeros(1,4);
for i = 1 : length(Ds)
   As{i} = pi*Ds{i}.^2/4;
   K_HP{i} = As{i}.^2 * 1e6/(8*pi*mu); % m^4 / MPa.s
   K_HP_sum(i) = sum(K_HP{i});
end
[gi_val, CN, CG, VAx] = gi(sim);

Rl = CN / mean(K_HP_sum); % MPa.s / m^4
[~,~,~,ktot,ktot_xa,~] =...
    compute_hydraulics(sim.gCond,1e6,0,...
    {sim.Dcross(1,:), sim.Dcross(2,:), sim.Dcross(3,:), sim.Dcross(4,:)},...
    VAx);

Rc = CN / ktot; %MPa.s / m^4

Rw = Rc - Rl;
Fw = Rw/Rc;