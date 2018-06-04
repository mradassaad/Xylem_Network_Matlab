# Xylem_Network_Matlab
This Matlab code is designed as a model to explore how water or emboli propagate throughout a flowering plant xylem tissue. This code is based on Object Oriented Programming (OOP) capabilities of Matlab. There are five objects that interact in this model:

* XylemNet: The xylem network object containing all conduits, conduit elements, ICCs, and clusters.

* Conduits: each conduit is made up of multiple consecutive conduit elements. Conduits are constrained to one direction, considered verticle along the orientation of water flow. One can think of conduits in this context as vessels.

* Conduit elements: condui elements form part of every conduit. All conduit elements are of the same length in the model. Therefore, the verticle distance between two consecutive nodes is the conduit element length which is set by the user.

* InterConduit Connections (ICCs): each ICC connects two conduits with each other. This is necessary because conduits do not usually span the whole length of the conductive tissue.

* Clusters: clusters are collections of conduits and ICCs that form an independent water pathway from the inlet of the conductive tissue to the outlet. A xylemNet object may contain multiple clusters.

There are two scripts a user may use to take advantage of the model without knowing its inner workings:

- [XylemNetIni.m](https://github.com/mradassaad/Xylem_Network_Matlab/blob/master/XylemNetIni.m) initializes a xylem network, pertinent parameters can be changed in that script. It allows you to edit the model parameters:
  - **rowNb**: the number of rows.
  - **colNb**: the number of columns.
  - **Pe**: probability that a node will terminate conduit element construction.
  - **NPe**: probability that a node will initiate conduit element construction.
  - **Pc**: probability that two horizontally adjacent nodes form an ICC.
  - **Lce**: conduit element length.
  - **Dc**: average conduit diameter.
  - **Dc_cv**: coefficient of variation of conduit diameters.
  - **Dp**: mean pit memebrane pore diameter.
  - **Dm**: mean pit membrane diameter.
  - **A**: scaling parameter of the Weibull distribution of pit membrane Bubble Propagation Pressures (BPPs).
  - **B**: curvature parameter of the Weibull distribution of pit membrane BPPs.
  - **fc**: contact fraction between two conduits.
  - **fpf**: pit field fraction.
  - **fap**: aperture area fraction.
  - **e_mean**: mean pit membrane strain at cavitation.
  - **e_cv**: coefficient of variation of pit membrane strain at cavitation.
  - **Tm**: pit membrane thickness.
  - **Lp**: pit chamber depth.
  - **BPPcalcmethod**: BPP calculation method (either <p style="color:DodgerBlue;">'Pore'</p> or <p style="color:DodgerBlue;">'Stretching'</p>)

- [VCGen.m](https://github.com/mradassaad/Xylem_Network_Matlab/blob/master/VCGen.m) allows you to compute the vulnerability curve associated with the generated xylem network.
