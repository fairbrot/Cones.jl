# Cones.jl

This package provides functionality related to polyhedral cones and polytopes.
It's main purpose is to provide utilties for the package [TailRiskScenGen.jl](https://github.com/STOR-i/TailRiskScenGen.jl).
It contains the following functionality:

- Provides type `PolyhedralCone` to represent cone defined by linear inequalities
- Provides type `FiniteCone` to represent cone defined by extremal rays
- Provides type `Polytope` to represent a polytope defined by extremal points and rays
- Implementation of Chernikova's algorithm which can be used to convert a `PolyhedralCone` to a `FiniteCone`
- Function to calculate extremal rays and points of a polytope defined by linear inequalities
- Functions to calculate projections of points onto polyhedral cones

# Algorithms

- N.V. Chernikova, Algorithm for finding a general formula for the non-negative solutions of a system of linear inequalities, USSR Computational Mathematics and Mathematical Physics, Volume 5, Issue 2, 1965, Pages 228-233, https://doi.org/10.1016/0041-5553(65)90045-5.
- Fernandez, Felipe & Quinton, Patrice. (1988). Extension of Chernikova's algorithm for solving general mixed linear programming problems. 

# Acknowledgements

This project makes uses code from the Siconos Numerics project
(http://siconos.gforge.inria.fr/) developed by INRIA for solving
linear complementarity complementarity problems using Lemke's
algorithm.