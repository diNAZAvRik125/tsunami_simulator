# Tsunami Simulator
## Physics Background:
The motion of fluid flow is described by the well-established Navier-Stokes equations. The shallow water equations are a simplifcation of the Navier-Stokes equations allowing for the motion of shallow layers of liquids under gravitational influence to be described. It should ne noted that the term "shallow" implies that the depth of the fluid is significantly smaller than its horizontal length. Some common applications of the shallow water equations include modelling tides and simulating tsunamis. 

The 1D shallow water equations are giving by:

![](https://github.com/a-rigido/Tsunami-Simulator/blob/master/README_images/SWEs.png)

where u is the filud velocity, n is the wave displacement, H is the depth of the floor, x is the horizontal space dimension (1D), and t is time. Converting these into Flux-Conservative form we get:

![](https://github.com/a-rigido/Tsunami-Simulator/blob/master/README_images/CF_form.png)

In Flux-Conservative form we can utilize the Two-Step Lax-Wendroff scheme (seen below) to solve for the vector U = [u,n] which gives us both the wave speed and height for the given simulation.

![](https://github.com/a-rigido/Tsunami-Simulator/blob/master/README_images/lax-wen.png)

So in essence, these scripts will utilize the 1D shallow water equations and the Two-Step Lax-Wendroff scheme to simulate swallow waves and a 1D tsunami.

## Highlights of Code:
* Utilizes Two-Step Lax-Wendroff Algorithm
* Visualization of 1D Shallow Water Equation for:
    * Topologically Flat Base (Unchanging Depth)
    ![](https://github.com/a-rigido/Tsunami-Simulator/blob/master/flat_bottom_SWE.png)        
    * Topologically Varying Base (Tsunami Simulation)
    ![](https://github.com/a-rigido/Tsunami-Simulator/blob/master/tsunami_simulation_SWE.png)
    

<p align="center">
Copyright (c) 2020 Alex Rigido
</p> 
