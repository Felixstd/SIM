# TODO List


1. Need to add initial conditions for mu-phi in ini_get. 
	Now, I seem to get infinities in the calculation which may be caused by that.
I need to output the pressure of the VP calculations.  

2. Periodic Boundary Conditions for the mu(I) - Phi(I) variables? 
	Need to look at what is done for pressure. 


3. Should put calculations near stress_strain


less 
The model converges more if using IMEX with the P = 0. Lets see how this goes. 
Work out with Ambrish the thing for the figures. 

IMEX seems to be not working at low mub but kinda converges at high mub. This is interesting.
I need to look at this with the 1d model. We could look at with IMEx and without. 
