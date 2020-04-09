DGV0D3V_direct_collision

Old copy of the code where convolution form of the collision operator is evaluated directly O(n^8) operations. Advantage is less aliasing and a beter use of resolution. This coud will be the course of training data for collision operator

Contributors: Alex Aleskeenko, Craig Euler

E-mail: alexander.alekseenko@csun.edu

This is a copy of a Fortran code that implements solution of the spatially homogeous Boltzmann equation using 
nodal-Discontinuous Galerkin discretization in velocity variable and Runge-Kutta/Adam-Bashforth integration in time. 
The code includes several models for evaluating the Boltzmann collision operator: the full Boltzmann collision operator 
for binary collisions (hard--spheres are currently supported), linearized collision operator, BGk, ES-BGK.

This is an older implementation of the method originally reported in 
2014 A.\ Alekseenko and E.\ Josyula, Deterministic solution of the spatially homogeneous Boltzmann equation using discontinuous Galerkin discretizations in the velocity space, \textit{Journal of Computational Physics}, Vol.~272, n.~1, (2014) 170--188.
2012 A.\ Alekseenko and E.\ Josyula, Deterministic solution of the Boltzmann equation using a discontinuous Galerkin velocity discretization.\ \textit{Proceedings of the $28^{\mathrm{th}}$ International Symposium on Rarefied Gas Dynamics, Spain 2012}, AIP Conference Proceedings, 2012, 8 pp.

The specifics of this code is that direct evaluaion of the 6D convolution form of the integral is used. Because the kernel in the convolution integral is sparse, the observed complexity of the algorithms is $n^8$, where $n$ is the number of points in one velocity dimension. The largest resolution tried for this code is 61^3 points. To make code tractable, it is MPI parallelized. 

Resently some important discoveires were made about the approach and a faster $O(n^6)$ algorithms are implemented and reported in the papers below:
2019 A. Alekseenko, J. Limbacher, Evaluating high order discontinuous Galerkin discretization of the Boltzmann collision integral in $O(N^2)$ operations using the discrete Fourier transform.\ \textit{Kinetic and Related Models} Vol.~12, n.~4,703-- 726, DOI: 10.3934/krm.2019027
2018 A.\ Alekseenko, T.\ Nguyen, and A.\ Wood, A deterministic-stochastic method for computing the Boltzmann collision Integral in $O(MN)$ operations. \textit{Kinetic and Related Models}, Vol., 11, no. 1937-5093

In the new papers Discrtete Fourier Transform is used to reduce the complexity of the formulation. However, this results to aliasing errors. Aliasing errors can be made small at the cost of using more grid points by padding solution with zeros. This effectively coasts resolution. 

This present code is the originally implemented version. It does not use FT and is not having aliasing errors. However, there is still a danger of truncating the domain too muich and other things that kinetic community is usually aware. 

While this code is impractical for solving Gas Dynamics Problems, it is practical for another use. We will use the code to produce trainign data for ANN approximations of collision operator. 

We will be developing a branch of the code, whose sole purpose will be to compute collision operatopr for a bunch of VDFs and save it on the hard drive. The computed outputs of the collision operator and the VDFs will become data for development of data driven formulatiosn fo the Boltzmann equation and for ANN approximations of collision operator. 
