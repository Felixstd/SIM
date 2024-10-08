ALL OF THE OUPUTS ARE IN output_sim_2

# Experiment 24
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 1
      Phi_0      = 1
      c_phi      = 1
with Phi_I is a constant = 1e-3


# Experiment 25
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 1
      Phi_0      = 1
      c_phi      = 1
with Phi_I is a constant = 0.5


# Experiment 26
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
IT's the base


# Experiment 27
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
phi = 0.5


# Experiment 28
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
phi = 0.05

# Experiment 29
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
With zeta capping being zetaC(i, j) = min(( mu_b + mu_I(i, j) / 2 ) * Pp(i, j) / shear_I(i, j), 2d8*Pstar)

# Experiment 30
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
With eta capping etaC(i, j)  = min((mu_I(i, j) / 2 ) * Pp(i, j) / shear_I(i, j), 1d10)



# Experiment 31
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
With zeta capping being zetaC(i, j) = min(( mu_b + mu_I(i, j) / 2 ) * Pp(i, j) / shear_I(i, j), 2d8)

# Experiment 32
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
With etaC(i, j)  = min((mu_I(i, j) / 2 ) * Pp(i, j) / shear_I(i, j), 1d15)



# Experiment 50
      d_average  = 1d3
      mu_0 = 3d-01
      mu_infty = 9d-01
      I_0        = 1e-3
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
      Without the use of IMEX. With a wind forcing of 10 m/s with tanh. 
      Now it seems to be converging. 


# Experiment 51
      d_average  = 1d4
      mu_0 = 3d-01
      mu_infty = 9d-01
      I_0        = 1e-3
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1
      Still converging, not much differences with before 

# Experiment 52
      d_average  = 1d3
      mu_0 = 3d-01
      mu_infty = 9d-01
      I_0        = 1e-3
      mu_b       = 10
      Phi_0      = 1
      c_phi      = 1

# Experiment 53
      d_average  = 1d3
      mu_0 = 3d-01
      mu_infty = 40d-01
      I_0        = 1e-3
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1


# Experiment 54
      d_average  = 1d3
      mu_0 = 3d-01
      mu_infty = 40d-01
      I_0        = 1e-1
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1

# Experiment 55
      d_average  = 1d3
      mu_0 = 3d-01
      mu_infty = 40d-01
      I_0        = 1e-1
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1

with only hibler p

# Experiment 56
      d_average  = 1d3
      mu_0 = 3d-01
      mu_infty = 40d-01
      I_0        = 1e-1
      mu_b       = 0.9
      Phi_0      = 1
      c_phi      = 1

with only hibler p






# CHANGED THE MU_PHI

The others were just important tests. 

The main thing to understand is that even with small time steps, I get convergence errors. 

# Experiment 94
      d_average  = 1d3*500
      mu_0 = 3d-01
      mu_infty = 6d-01
      I_0        = 0.3
      mu_b       = 1
      Phi_0      = 1
      c_phi      = 1
      delta_x = 1km

# Experiment 96
      d_average  = 1d3
      mu_0 = 1d-01
      mu_infty = 9d-01
      I_0        = 1d-03   
      mu_b       = 100
      Phi_0      = 1
      c_phi      = 1
With hibler pressure