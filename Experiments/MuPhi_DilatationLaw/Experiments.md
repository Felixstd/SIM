# Experiment 1



This is the control experiment made with the mu(I)-phi(I)
rheology without the dilatancy law.

This was at Nx, Ny = 400, 1000 with delta_x = 2km and delta_t = 0.2s. 

This is also done with the Picard solver wihtout the use of IMEX. 

With a convergence criterion of 1d-15 and 2000 maximum loops. 

With the following parameters:
d_average  = 1d3
I_0        = 1e-3
mu_0 =     0.1
mu_infty = 0.9
mu_b       = 1
Phi_0      = 1
c_phi      = 1

I think I could use a smaller delta_x, but it's good enough for now. 

## Problems with 1
Initially, I had set that for h < 1e-6, the friction coefficient to be mu_0. 
But, I think it should be mu_infty rather. Because open water refers to I -> infinity. 
I will rerun it and see what happens. 



# Experiment 2
Same as Experiment 1 but with the new initial conditions. 


There seems to be no real differences. 

# Experiment 3

This was made with the mu(I)-Phi(I) rheology with the dilatancy law. 

This is at Nx, Ny = 400, 1000 with delta_x = 2km and delta_t = 0.2s.

With a convergence criterion of 1d-15 and 2000 maximum loops.


Also, I used the upwind integration for the continuity equations. 

With the following parameters:
d_average  = 1d3
I_0        = 1e-3
mu_0 =     0.1
mu_infty = 0.9
mu_b	   = 1
Phi_0	   = 1
c_phi	   = 1
K = 2


This simulation as is is giving not what we expect. But, It's not exploding even after 
30 seconds of simulation. We thus observe that for now, the velocities are too high and that the 
ice is not breaking. 


For this, I used the same $\eta$ has for the mu(I)-Phi(I) case and $\zeta = 0$, but I don't 
think it's the way to go. I wanted to isolate both components. 


The ice does not seem to deform at the same time, the ice starts deforming at the top but as time passes
the deformation king of moves towards the lower boundary. But, there's no breaking. 


## Problems 
From what I see, it should be more something like 
$\tau = \mu_{eff} p_{eq} = \mu(I) + \tan \psi$ with 
$\tan \psi = K(A - phi)$. 

So here, $ eta = \mu_{eff} p_{eq}$ and $\zeta = 2\eta$.


# Experiment 4

Same parameters as 2 but with 
$\zeta = 2\eta = 2\mu_{eff}(I)p$. 
And the new initial conditions. 

Here, again we see that as time goes, the deformation extends and wants to go to the bottom of the ice. 
But, interestingly here, a triangle at the top starts to appearing. The deformation is full of holes, 
which is something I need to look at, but we have a nice triangle at the top and with which seems a ~50° angle. 


Also, I did not regularize the viscous coefficients here. Because P is already.

But the alignment condition says that 
$\frac{\tau}{\norm{\tau}} = \frac{D}{\norm{D}}$. 

And, here I don't have the shear. 

# Experiment 5

Same parameters as 2 but now with


$\zeta = 2\eta = 2\mu_{eff}(I)p/eII$.

Here, it seems that we obtain some type of fracture features as time goes on. 
But, again, it's weird that the ice does not seem to break at the bottom and that it takes
time for the deformation to go all the way down. Maybe the parameters? 

Or the way we are computing the stresses and divergence. 


# Experiment 6
Let's do the same thing but with the fact that eta = zeta = mu(I)p/eII. 
So, this corresponds to the mu(I) rheology. 


This is done with the semilag advection scheme. 

ATTENTION, I forgot the $\tan \psi$ term in the mu(I) calculation. 

I need to add that. I think it looks kind of the same as the other ones because of that 
especially.  

Interesting!!

This seems to work better than all the others.

It looks a lot like 2!! 
Actually, it seems to be the same but with less places where we have weird stuff happening

Would need to do a sensitivity analysis to see how this is impacted by the different parameters.

# Experiment 7

Same as 6 but with the right mu. 


Weird that I don't have any thing for A and h. It doesnt seem to evolve. 


Because because of the advection scheme? 

# Experiment 8
Same as 7 but with an upwind advection scheme. 


# Notes

In all of these abose, the pressure used is the regularized one from mu_phi. I would like to test if 
I only use the one in Hibler.


# Expriment 9
Same as all of the other but rather using p_eq I use P_max only.




## Sensitivity analysis
	
