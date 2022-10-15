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

## Overview of the solvers
&emsp;There are totally 6 solvers in **CP3d**: 

<img src="./doc/Overview-6-solvers.png" width="40%" height="40%" div align=center />


## Acknowledgements :clap:
&emsp;Since Sep 2019, when I finally decided to develop my own CFD-DEM code from scratch, I have learnt quite a lot from the following really kind researchers (**in alphabetical sequence**):

* [Dr. Costa](https://p-costa.github.io/) from University of Iceland, and his second-order DNS code [CaNS](https://github.com/p-costa/CaNS), also his papers on IBM approach.
* [Dr. He](https://www.engineering.iastate.edu/people/profile/phe/) from Iowa State University, and his fourth-order DNS solver [HercuLES](https://github.com/friedenhe/hercules).
* [Prof. Ji](http://faculty.tju.edu.cn/ChunningJi/en/index.htm) from Tianjin University, on the fruitful discussion about the particle IBM method, and on the access to their in-house DNS/LES-Solid interaction code **_cgLES_**.
* [Dr. Laizet](http://www.imperial.ac.uk/people/s.laizet) from Imperial College London, and their compact FD code [Incompact3d](https://github.com/xcompact3d/Incompact3d).
* [Prof. Marchioli](http://158.110.32.35/) from University of Udine, on the fruitful and continuous discussion about one-way CFD-Particle coupling benchmark and on the access to their [benckmark data](http://158.110.32.35/download/DNS-TEST-CASE/).
* [Prof. Meiburg](https://me.ucsb.edu/people/eckart-meiburg) from University of California, Santa Barbara.
* [Dr. Norouzi](https://www.researchgate.net/profile/Hamid-Norourzi) from University of Tehran, and his book **_Coupled CFD‐DEM Modeling: Formulation, Implementation and Applimation to Multiphase Flows_**, besides the [attached DEM code](https://www.wiley.com//legacy/wileychi/norouzi/form.html?type=SupplementaryMaterial).
* [Prof. Orlandi](http://dma.ing.uniroma1.it/users/orlandi/resume.html) from Sapienza University of Rome, and his book **_Fluid flow phenomena: a numerical toolkit_**, besides the [attached CFD code](http://dma.ing.uniroma1.it/users/orlandi/diskette.tar.gz).
* [Dr. Tschisgale](https://www.researchgate.net/profile/Silvio-Tschisgale) from Institute of Air Handling and Refrigeration, on the fruitful and continuous discussion about their IBM approach.
* [Prof. Zhao](http://www.hy.tsinghua.edu.cn/info/1154/1829.htm) from Tsinghua university, on the one-way CFD-Particle coupling benchmark.
* ......

## Contact and Feedback :email:
&emsp;If you have any question, or want to contribute to the code, please don't hesitate to contact me: Zheng Gong (gongzheng_justin@outlook.com)
