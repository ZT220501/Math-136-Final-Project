# Math-136-Final-Project
This project is a simulation of the unidirectional flow problem and flow in a cavity problem.

![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/uniflow_explicit_trivial_boundary.gif) 
![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/uniflow_explicit_trivial_boundary_unstable.gif)

The two animations above are the explicit finite difference results for unidirectional flow horizontal velocity as a function of vertical axis $y$. The first one is the caes when $\alpha=0.49$ and the second one is the case when $\alpha=0.51$. Clearly the case $\alpha=0.51$ is extremely unstable.

![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/uniflow_implicit_trivial_boundary0.49.gif)
![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/uniflow_implicit_trivial_boundary0.51.gif)

The two animations above are the implicit finite difference results for unidirectional flow horizontal velocity as a function of vertical axis $y$. The first one is the caes when $\alpha=0.49$ and the second one is the case when $\alpha=0.51$. In both case, the finite difference method gives stable result.

Results of unidirectional flow non-trivial constant boundary condition and sinuosoidal boundary conditionss can be found under the Result folder.

![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/Low_Reynold_Simulation_Re2.gif)
![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/Middle_Reynold_Simulation_Re200.gif)
![](https://github.com/ZT220501/Math-136-Final-Project/blob/main/Result/High_Reynold_Simulation_Re2000.gif)

The three animations above are the simulations of the fluid velocity field in the order when Reynolds number $Re=2$, $Re=200$, and $Re=2000$. The steady cases can be found in the Result folder.
