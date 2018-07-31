%rowNb and colNb are respectively the number of rows and columns of the Xylem Network (XN).
rowNb = 100;
colNb = 1000;
%Pc and NPc are respectively the probability of initiating and terminating a conduit.
%Once a conduit is initiatied, the model will keep constructing consecutive conduit elements
%until a terminating node is reached with probability NPc.
%Pc and NPc are tuned so that conduit lengths and connectivity match real measurements.
Pc = 0.91;
NPc = 0.87;
%Pe is the probability of two horizontally adjacent nodes, forming parts of two
%different conduits, will form a InterConduit Connection (ICC).
%Pe is tuned so that conduit connectivity matches real measurements.
Pe = 0.5;
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
%BPPcalcmethod below is set to 'Stretching'.
e_mean = 0.02;
e_cv = 0.35;
%Tm is the average thickness of membranes in meters.
Tm = 234e-9;
%Lp is the average pit chamber depth in meters.
Lp = 654e-9;
%BPPcalcmethod specifies whether the Bubble Propagation Pressure(BPP) is to be 
%prescribed by either the width of the largest pore or  membrane stretching.
%Should be set to either 'Pore' or 'Stretching'.
BPPcalcmethod = 'Pore';

%sim will be the XylemNet object created using the parameters defined above
sim = XylemNet(rowNb,colNb,Pc,NPc,Pe,Lce,Dc,Dc_cv,Dp,Dm,A,B,fc,fpf,fap,...
      e_mean,e_cv,Tm,Lp,BPPcalcmethod);
