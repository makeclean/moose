# RBE3Constraint

!syntax description /Constraints/RBE3Constraint

The RBE3 (Rigid Body Element 3) constraint constrains slave nodes to move as a weighted combination of master nodes. This is commonly used in structural mechanics to distribute loads or motions from a reference point to a set of nodes on a boundary.

The constraint supports three weighting types:
- `constant`: Equal weights are assigned to all master nodes regardless of slave position
- `inverse_distance`: Weights are computed based on inverse distance from each slave node to each master node, giving more influence to closer masters
- `manual`: User-specified weights for each master node

!syntax parameters /Constraints/RBE3Constraint

!syntax inputs /Constraints/RBE3Constraint

!syntax children /Constraints/RBE3Constraint

!bibtex bibliography
