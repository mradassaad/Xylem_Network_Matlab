# Xylem_Network_Matlab
This Matlab code is designed as a model to explore how water or emboli propagate throughout a flowering plant xylem tissue. This code is based on Object Oriented Programming (OOP) capabilities of Matlab. There are five objects that interact in this model:
1- XylemNet: The xylem network object containing all vessels, vessel elements, IVCs, and clusters.
2- Vessels: each vessel is made up of multiple consecutive vessel elements. Vessels are constrained to one direction, considered verticle along the orientation of water flow.
3- Vessel elements: vessel elements form part of every vessel. All vessel elements are of the same length in the model. Therefore, the verticle distance between two consecutive nodes is the vessel element length which is set by the user.
4- InterVessel Connections (IVCs): each IVC connects two vessels with each other. This is necessary because vessels do not usually span the whole length of the conductive tissue.
5- Clusters: clusters are collections of vessels and IVCs that form an independent water pathway from the inlet of the conductive tissue to the outlet. A xylemNet object may contain multiple clusters.

There are two scripts a user may use to take advantage of the model without knowing its inner workings:
XylemNetIni.m initializes a xylem network, pertinent parameters can be changed in that script. It allows you to edit the model parameters.
