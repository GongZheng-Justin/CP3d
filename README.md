# CP3d
&emsp;**CP3d** (Channel-Particle 3d) is a comprehensive Euler-Lagrange solver for the direct numerical simulations of particle-laden flows.

* Have the ability to perform **one-way**, **two-way**, **interface-unresolved four-way**, and **interface-resovled four-way** coupling regime.
* For the fluid sub-solver, both the second-order and four-order finite difference discretization are available, in company with three viscous term treatment approaches: i.e., fully-implicit, partial-implicit, and fully-explicit.
* For the **discrete-element-method (DEM)** solver, both linear and non-linear contact force models can be used, and the collision time can be also adaptive.
* The **immersed boundary method (IBM)** is used in interface-resolved simulation, and totally three IBM coupling approaches are included.
* In order to improve the numerical accuracy of the computational **Basset history force**, a third-order exponential approximation method is proposed in interface-unresolved four-way regime.
* An **averaged lubrication force model** is proposed for the short-range hydrodynamic interaction.
* The volume integration approach is also modified to adapt the staggered mesh configuration.
* The resulting solver is able to simulate large scale cases of **billions of grid points** with **millions of moving particles** in interface-resolved four-way regime using only hundreds of computational cores.

## Contact and Feedback :email:
&emsp;If you have any question, or want to contribute to the code, please don't hesitate to contact me: Zheng Gong (gongzheng_justin@outlook.com)
