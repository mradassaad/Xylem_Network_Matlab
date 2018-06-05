# Xylem_Network_Matlab

## Introduction
Research on plant hydraulics is advancing at great pace and from different scientific angles. There is a need to combine newly acquired knowledge on plant hydraulic tissue anatomy, physiology, and whole tissue hydraulic function. This model was developped to do just that for flowering plant xylem. The object oriented nature of the model provides modularity and allows for easy implementation of additional levels of complexity of xylem function.

The Xylem Network model first builds a xylem network (XN) of random topology. Topology means the relative position of vessels relative to each other. Properties such as mean vessel length, diameter and mean mean pit membrane size are controlled through input parameters to the model as detailed in the [How To Use section](#how-to-use).

The model can then generate Vulnerability Curves (VCs) and extract hydraulic properties of the generated XN. Below is a figure illustrating how emboli propagate through the xylem as the pressure difference between the water and the atmosphere (![alt text][DP]) increases.

![alt text][cav_gif]

## Classes
This Matlab code is designed as a model to explore how water or emboli propagate throughout a flowering plant xylem tissue. This code is based on Object Oriented Programming (OOP) capabilities of Matlab. There are five objects that interact in this model:

* XylemNet: The xylem network object containing all conduits, conduit elements, ICCs, and clusters.

* Conduits: each conduit is made up of multiple consecutive conduit elements. Conduits are constrained to one direction, considered verticle along the orientation of water flow. One can think of conduits in this context as vessels.

* Conduit elements: condui elements form part of every conduit. All conduit elements are of the same length in the model. Therefore, the verticle distance between two consecutive nodes is the conduit element length which is set by the user.

* InterConduit Connections (ICCs): each ICC connects two conduits with each other. This is necessary because conduits do not usually span the whole length of the conductive tissue.

* Clusters: clusters are collections of conduits and ICCs that form an independent water pathway from the inlet of the conductive tissue to the outlet. A xylemNet object may contain multiple clusters.

## How to use

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
  - **BPPcalcmethod**: BPP calculation method (either <font color="blue">'Pore'</font> or <font color="blue">Stretching'</font>)

- [VCGen.m](https://github.com/mradassaad/Xylem_Network_Matlab/blob/master/VCGen.m) allows you to compute the vulnerability curve associated with the generated xylem network.


[DP]: http://www.sciweavers.org/tex2img.php?eq=%5CDelta%20P&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0
[cav_gif]: https://github.com/mradassaad/Xylem_Network_Matlab/blob/master/Images/cav_gif.gif
