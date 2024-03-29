%rowNb and colNb are respectively the number of rows and columns of the Xylem Network (XN).
rowNb = 50;
colNb = 10;
depNb = 50;
%Pc and NPc are respectively the probability of initiating and terminating a conduit.
%Once a conduit is initiatied, the model will keep constructing consecutive conduit elements
%until a terminating node is reached with probability NPc.
%Pc and NPc are tuned so that conduit lengths and connectivity match real measurements.
Pc = [0.75 0.75];
NPc = [0.75 0.75];
%Pe is the probability of two horizontally adjacent nodes, forming parts of two
%different conduits, will form a InterConduit Connection (ICC).
%Pe is tuned so that conduit connectivity matches real measurements.
Pe_rad = [0.9 0.9];
Pe_tan = [0.02 0.02];
%radDist: It interpolates Pc and Pe values between vessel closest to the pith
%(index 1) and those farthest (last element).
% radDist = [1 1 1  0.1 0.1 0.1 0.1 0.1 0.1 0.01];
radDist = ones(1, colNb);
%Lce is the length of a single conduit elements in meters
Lce = 2.88e-3;
%Dc is the average conduit diameter in meters
Dc = 24.6e-6;
%Dc_cv is the coefficient of variation of conduit diameters
Dc_cv = 0.12;
%Dp is the average pore diameter in meters
Dp = 32.0e-9;
%Dm is the average membrane diameter in meters
Dm = 6.3e-6;
%A and B are the parameters of the Weibull CDF from which large pore diameters
%are sampled. This will only be use if BPPcalcmethod below is set to 'Pore'.
A = 16;
B = 4;
%fc and fpf are the average contact fraction and pit field fraction between 
%two conduits, respectively.
fc = 0.31;
fpf = 0.7;
%fap is the aperture fraction.
fap = 0.06;
%e_mean and e_cv are the mean and coefficient of variation of the normal distribution
%from which membrane stretching limits are sampled. These are only used when 
%ASPcalcmethod below is set to 'Stretching'.
e_mean = 0.02;
e_cv = 0.35;
%Tm is the average thickness of membranes in meters ('Stretching').
Tm = 234e-9;
%Lp is the average pit chamber depth in meters ('Stretching').
Lp = 654e-9;
%ASPcalcmethod specifies whether the Air Seeding Pressure(ASP) is to be 
%prescribed by either the width of the largest pore or  membrane stretching.
%Should be set to either 'Pore' or 'Stretching'.
ASPcalcmethod = 'Pore';
%Note: the 'Pore' ASPcalcmethod can also be used as a random ASP pore
%generator to match VCs to each other.

%sim will be the XylemNet object created using the parameters defined above
sim = XylemNet(rowNb,colNb,depNb,Pc,NPc,Pe_rad,Pe_tan,radDist, Lce,Dc,Dc_cv,Dp,Dm,A,B,fc,fpf,fap,...
      e_mean,e_cv,Tm,Lp,ASPcalcmethod);
