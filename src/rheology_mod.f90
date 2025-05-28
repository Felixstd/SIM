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


MODULE muphi

! mu_0       : 
! mu_infty   : 
! c_phi      : 
! I_0        : 
! d_average  : 
! phi_0      :

  IMPLICIT NONE

  LOGICAL :: devstrain, Water_Col, A2Phi
  LOGICAL :: step_water, Phi_eq
  DOUBLE PRECISION :: mu_0, mu_infty 
  DOUBLE PRECISION :: c_phi
  DOUBLE PRECISION :: I_0
  DOUBLE PRECISION :: d_average
  DOUBLE PRECISION :: Phi_0, Phi_max, Phi_min
  DOUBLE PRECISION :: tau
  DOUBLE PRECISION :: mu_b
  DOUBLE PRECISION :: K_div
  DOUBLE PRECISION :: c_1, c_2
  DOUBLE PRECISION :: phi_f_micro
  DOUBLE PRECISION :: theta



END MODULE muphi

