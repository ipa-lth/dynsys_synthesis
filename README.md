# Automatic Synthesis of a Saturating State Space Controller Based on Convex Optimization for Task Space Controlled Industrial Robots
#### Lorenz Halt , Philipp Tenbrock , Frank Nägele and Andreas Pott 

## Software contribution to the submitted research paper

This repository holds all necessary data, including the identified model in form of a npy-file (numpy pickle) and the convex optimisation problem itself.

The [jupyter notebook](Dynamical_system_control_synthesis.ipynb) loads all data and executes an optimisation. To run it on your own computer please make sure to resolve all dependencies (python 2.7):
* numpy
* cvxpy (incl. solver: CVXOPT, MOSEK(1), SCS)
* matplotlib





(1) A licence is necessary for MOSEK. However, any suitable SDP solver can be included.
