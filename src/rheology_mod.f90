MODULE ellipse
!
! C : ice strength parameter
! Pstar : Ice compression strength parameter
! ell2  : ellipticity**2 
! ell_2 : 1/ellipticity**2
  IMPLICIT NONE
  DOUBLE PRECISION :: C, Pstar, ell2, ell_2, Tens

END MODULE ellipse



MODULE triangle
!
! Cohe : Ice tensile strength
! phi : internal angle of friction
! delta : angle of dilatancy
! etamax : maximum shear viscosity
  IMPLICIT NONE
  DOUBLE PRECISION :: Cohe, phi, delta, etamax

END MODULE triangle



MODULE elastic

! phi : internal angle of friction
! sigC : tensile strength cut-off
! Poisson : Poisson ratio of sea ice
! Young : Young's ration of sea ice 
! lambda0: viscous relaxation time scale

  IMPLICIT NONE

  DOUBLE PRECISION :: Young, Poisson, lambda0
  DOUBLE PRECISION :: alpha, Theal 
  DOUBLE PRECISION :: Cohe, sigt, sigc, phi, Tdam
  DOUBLE PRECISION :: theta_cor

END MODULE elastic

MODULE MC_teardrop

! Parameters for the MC teardrop rheology
! mu_0       : quasi-static friction coefficient 
! mu_inf     : kinetic friction coefficient 
! I0         : transition inertial number 
! gamma_kt   : max of yield curve 
! alpha_reg  : regularization constant 

  IMPLICIT NONE 

  DOUBLE PRECISION :: mu_0, mu_inf, I0
  DOUBLE PRECISION :: gamma_kt, cyc
  DOUBLE PRECISION :: alpha_reg
  DOUBLE PRECISION :: phi_prime
  DOUBLE PRECISION :: delta_c
  DOUBLE PRECISION :: zetamax
  DOUBLE PRECISION :: repl_press
  INTEGER :: non_normal

END MODULE MC_teardrop
