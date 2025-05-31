! This collection of subroutines is used for the 
! mu-phi rheology. It contains all of the necessary 
! subroutines to compute phi, I, shear, mu, psi. 
!
! Written by Félix St-Denis, 2024-2025


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing different                         !
!                  quantities at one time step.                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inertial_number

!------------------------------------------------!
!--- SUBROUTINE COMPUTING THE INERTIAL NUMBER ---!
!------------------------------------------------!

    !-- Modules to include --!
    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    !-- Local variables --!

    integer i, j, peri ! indexes
    double precision highI, lowI, Press ! Max I, min I, Pressure(i, j)

    ! Need to look for the eps and see if it is the right thing to do. 
    ! Also, for the value of the High I as P -> 0

    !-- Setting values --!
    ! highI = 1d0
    highI = 1d0
    lowI  = 0d0

    !-- Loop through all the grid cells --!
    do i = 0, nx+1
        do j = 0, ny+1

            ! if ocean 
            if (maskC(i,j) .eq. 1) then

                    ! Capping the pressure to eliminate singularities
                    Press = max(Pp(i, j), 1d-20)
                    if (devstrain) then 
                        inertial(i, j) = min(SQRT(rhoice * h(i, j)/Pp(i, j)) * d_average*2d0*shearC_I(i, j), highI)
                    else
                        ! Computing inertial number
                        if (h(i, j) < 1e-3) then    
                            inertial(i,j) = 1d0
                        else
                            inertial(i, j) = min(SQRT(rhoice * h(i, j)/Pp(i, j)) * d_average*shearC_I(i, j), highI)
                        endif
                    endif
            endif
        enddo
    enddo

    !-- Boundary Conditions --!


    if (Periodic_y .eq. 0) then
        do i = 0, nx+1

            if (maskC(i,0) .eq. 1) then             
                inertial(i,1)  = lowI
            endif

            if (maskC(i,ny+1) .eq. 1) then             
                inertial(i,ny)  = lowI
            endif

        enddo
    endif            
	  
    if (Periodic_x .eq. 0) then
        do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then             
                inertial(1,j)  = lowI
            endif

            if (maskC(nx+1,j) .eq. 1) then             
                inertial(nx,j)   = lowI  
            endif           

        enddo
    endif
    
end subroutine inertial_number

subroutine angle_friction_mu 

    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j, peri
    double precision eps, scaled_A, diff_A, min_inertial

    peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions      

    eps = 1d-12

    do i = 0, nx+1
        do j = 0, ny+1

            if (maskC(i,j) .eq. 1) then
                    mu_I(i, j) = mu_0 +  ( mu_infty - mu_0 ) / ( (I_0*c_phi)/inertial(i, j) + 1 )

            endif
        
        enddo
    enddo

    if (Periodic_y .eq. 0) then
        do i = 0, nx+1

            if (maskC(i,0) .eq. 1) then             
                mu_I(i,1)  = mu_0
            endif

            if (maskC(i,ny+1) .eq. 1) then             
                mu_I(i,ny)  = mu_0
            endif

        enddo
    endif            
	  
    if (Periodic_x .eq. 0) then
        do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then             
                mu_I(1,j)  = mu_0
            endif

            if (maskC(nx+1,j) .eq. 1) then             
                mu_I(nx,j)   = mu_0  
            endif           

        enddo
    endif

    

end subroutine angle_friction_mu

subroutine shear_inv(utp, vtp)

    use datetime, only: datetime_type
    use ellipse
    use muphi
    
    implicit none

    include 'parameter.h'
    include 'CB_Dyndim.h'
    include 'CB_DynVariables.h'
    include 'CB_DynForcing.h'
    include 'CB_const.h'
    include 'CB_mask.h'
    include 'CB_options.h'

    integer i, j, peri, ii, jj, im, jm
    character filename*64

    double precision dudx, dvdy, dudy, dvdx, land, lowA, summaskC
    double precision dudy1(0:nx+2,0:ny+2), dudx1(0:nx+2,0:ny+2), dvdy1(0:nx+2,0:ny+2), dvdx1(0:nx+2,0:ny+2)
    double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

    peri = Periodic_x + Periodic_y

    land  = -999d0
    ! land  = 0d0
    shearC_I = land
    lowA = 0d0

    if (peri .ne. 0) call periodicBC(utp,vtp)

    dudx       = 0d0
    dvdy       = 0d0
    dudy       = 0d0
    dvdx       = 0d0
    
    !---- At the grid centers (C-Grid) ----!
    !--- Only for the domain ---!
    do i = 1, nx
        do j = 1, ny
            
            if ( maskC(i,j) .eq. 1 ) then
                    
                  
                    dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                    dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax

                  
                if (( maskC(i+1,j) + maskC(i-1,j) .eq. 2 )) then
                    dvdx = ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                            ( vtp(i-1,j) + vtp(i-1,j+1) ) ) /      &
                            ( 4d0 * Deltax )
                     
                elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. 1 ) then
                    dvdx = ( 1d0 * ( vtp(i+1,j) + vtp(i+1,j+1) ) +  &
                            3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) /  &
                            ( 6d0 * Deltax )
                
                elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. -1 ) then
                    dvdx = ( -1d0 * ( vtp(i-1,j) + vtp(i-1,j+1) ) - &
                            3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) / &
                            ( 6d0 * Deltax )
                    
                elseif ( maskC(i+1,j) + maskC(i-1,j) .eq. 0 ) then
                    
                    print *, 'WARNING: irregular grid cell case1', i, j
                    
                endif
               
                if     ( maskC(i,j+1) + maskC(i,j-1) .eq. 2 ) then
                    !---- For open boundary condition at the top ----!
                    if (j .eq. ny) then 

                        if (Periodic_y .eq. 0) then
                            dudy = (utp(i, j) + utp(i+1, j) - &
                                   (utp(i, j-1) + utp(i+1, j-1)))/(4d0*Deltax)
                    
                        endif
                    else
                        dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                            ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                            ( 4d0 * Deltax )
                    endif
                     
                elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. 1 ) then
                    dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                            3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                            ( 6d0 * Deltax )
                     
                elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. -1 ) then
                    dudy = ( -1d0 * ( utp(i,j-1) + utp(i+1,j-1) ) - &
                            3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                            ( 6d0 * Deltax )
                     
                elseif ( maskC(i,j+1) + maskC(i,j-1) .eq. 0 ) then
                     
                    print *, 'WARNING: irregular grid cell case2',i,j
                     
                endif

!----- stresses and strain rates at the grid center -------------------------   

                if (devstrain) then

                    shearC_I(i,j) = sqrt(1/2d0*(dudx**2d0 + dvdy**2d0 + &
                                    2d0*(1/2d0*(dudy+dvdx)**2) - 1/2d0*(dudx + dudy)**2d0))
                    if ((dudx**2d0 + dvdy**2d0 + &
                                    2d0*(1/2d0*(dudy+dvdx)**2)) < ((dudx + dudy)**2d0)) then
                        print*, 'shear problem'
                    endif
                else 
                    shearC_I(i,j) = sqrt((dudx - dvdy )**2d0 &  
                            + ( dudy + dvdx )**2d0)
                endif
            else
                shearC_I(i,j) = land

            endif

        enddo
    enddo
    

    !------ Boundary Conditions ------!
    ! print*, 'MASKC' , maskC(0,ny+1)
    !----- Boundary Conditions in X ---!
    if (Periodic_x .eq. 1) then
        
        if (Periodic_y .eq. 0) then
            
            !------- i = 0 --------!
            i = 0
            do j = 1, ny

                if (maskC(i,j) .eq. 1) then
                    dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                    dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax

                    dvdx =  ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                            ( vtp(nx-1,j) + vtp(nx-1,j+1) ) ) /      &
                            ( 4d0 * Deltax )

                    !----- dudy (i = 0)------!
                    if (j.eq. 1) then 
                    !------ This is assuming a wall at the bottom -------!
                        dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                            3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                            ( 6d0 * Deltax )
                    
                    else 
                        if ( maskC(i,j+1) + maskC(i,j-1) .eq. 2 ) then
                            !---- For open boundary condition at the top ----!
                            if (j .eq. ny) then 

                                if (Periodic_y .eq. 0) then
                                    dudy = (utp(i, j) + utp(i+1, j) - &
                                            (utp(i, j-1) + utp(i+1, j-1)))/(4d0*Deltax)
                    
                                endif
                            else
                                dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                                        ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                                        ( 4d0 * Deltax )
                            endif

                        elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. 1 ) then
                            dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                                3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                                ( 6d0 * Deltax )
                        
                        elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. -1 ) then
                            dudy = ( -1d0 * ( utp(i,j-1) + utp(i+1,j-1) ) - &
                                3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                                ( 6d0 * Deltax )
                        endif
                    endif 

                    if (devstrain) then

                        shearC_I(i,j) = sqrt(1/2d0*(dudx**2d0 + dvdy**2d0 + &
                                        2d0*(1/2d0*(dudy+dvdx)**2) - 1/2d0*(dudx + dudy)**2d0))
                        if ((dudx**2d0 + dvdy**2d0 + &
                                        2d0*(1/2d0*(dudy+dvdx)**2)) < ((dudx + dudy)**2d0)) then
                            print*, 'shear problem'
                        endif
                    else 
                            shearC_I(i,j) = sqrt((dudx - dvdy )**2d0 &  
                                + ( dudy + dvdx )**2d0)
                    endif
                else

                    shearC_I(i,j) = land

                endif
            enddo

            !-------- i = nx+1 --------!
            i = nx+1
            do j = 1, ny

                if (maskC(i,j) .eq. 1) then
                    dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                    dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax

                    dvdx =  ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                            ( vtp(i-1,j) + vtp(i-1,j+1) ) ) /      &
                            ( 4d0 * Deltax )

                    !----- dudy -----!

                    if (j .eq. 1) then
                        
                        !----- Forward Spatial Derivative at the Boundary -----!
                        
                        dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                            3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                            ( 6d0 * Deltax )

                    else 


                        if (( maskC(i,j+1) + maskC(i,j-1) .eq. 2 )) then
                            !---- For open boundary condition at the top ----!
                            if (j .eq. ny) then 

                                if (Periodic_y .eq. 0) then
                                    dudy = (utp(i, j) + utp(i+1, j) - &
                                            (utp(i, j-1) + utp(i+1, j-1)))/(4*Deltax)
                    
                                endif
                            else
                                dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                                        ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                                        ( 4d0 * Deltax )
                            endif


                        elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. 1 ) then
                            dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                                3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                                ( 6d0 * Deltax )
                        
                        elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. -1 ) then
                            dudy = ( -1d0 * ( utp(i,j-1) + utp(i+1,j-1) ) - &
                                3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                                ( 6d0 * Deltax )
                        endif

                    endif

                    if (devstrain) then

                        shearC_I(i,j) = sqrt(1/2d0*(dudx**2d0 + dvdy**2d0 + &
                                        2d0*(1/2d0*(dudy+dvdx)**2) - 1/2d0*(dudx + dudy)**2d0))
                        if ((dudx**2d0 + dvdy**2d0 + &
                                        2d0*(1/2d0*(dudy+dvdx)**2)) < ((dudx + dudy)**2d0)) then
                            print*, 'shear problem'
                        endif
                    else 
                            shearC_I(i,j) = sqrt((dudx - dvdy )**2d0 &  
                                + ( dudy + dvdx )**2d0)
                    endif
                else
                    shearC_I(i,j) = land
                
                endif

            enddo
        endif

    else 
        !--- Open Boundary Conditions (Neumann BC) in x ---!
        i = 0
        do j = 1, ny 

            if (maskC(i,j) .eq. 1) then
                dudx = 0d0
                dvdy = 0d0 
                dvdx = 0d0
                dudy = 0d0
            else
                shearC_I(i,j) = land
            endif
        enddo

        i = nx+1
        do j = 1, ny 

            if (maskC(i,j) .eq. 1) then
                dudx = 0d0
                dvdx = 0d0
                dvdy = 0d0
                dudy = 0d0

            else
                shearC_I(i,j) = land
            endif
        enddo


    endif

    !--- Boundary Conditions in Y ---!
    
    if (Periodic_y .eq. 1) then
    !--- Periodic (Dirichlet BC) in y ---!
    !------- Not debugged -------!
        print*, 'Stop not implemented'
        stop

    endif

    return
end subroutine shear_inv


subroutine volumefraction_phi(An_1, Aout)

!----------------------------------------------------------!
!--- SUBROUTINE USED TO COMPUTE THE VOLUME FRACTION PHI ---!  
!---   FROM THE EMPIRICAL LAW DERIVED FROM EXPERIMENTS  ---!    
!----------------------------------------------------------!

!--- Needs to be run after calling inertial_I() ---!
    
    !-- Loading modules --!
    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    !-- Local Variables --!
    integer i, j, peri
    double precision eps
    double precision                :: An_1(0:nx+2,0:ny+2)
    double precision, intent(out)   :: Aout(0:nx+2,0:ny+2)

    eps = 1d-08
    
    if (peri .ne. 0) then
        call periodicBC2(An_1)
    endif

    if (Periodic_y .eq. 0) then

        do i = 0, nx+1
               
            if (maskC(i,0) .eq. 1) then

                An_1(i,0) = ( 4d0 * An_1(i,1) - An_1(i,2) )/3d0
                An_1(i,0) = max(An_1(i,0), 0d0)
                An_1(i,0) = min(An_1(i,0), 1d0)
                  
             endif

             if (maskC(i,ny+1) .eq. 1) then

                An_1(i,ny+1)= ( 4d0 * An_1(i,ny) - An_1(i,ny-1) ) / 3d0
                An_1(i,ny+1)= max(An_1(i,ny+1), 0d0)
                An_1(i,ny+1)= min(An_1(i,ny+1), 1d0)

             endif
 
          enddo
                  
      endif            
	  
      if (Periodic_x .eq. 0) then
	  
         do j = 0, ny+1

             if (maskC(0,j) .eq. 1) then
                 An_1(0,j)  = ( 4d0 * An_1(1,j) - An_1(2,j) ) / 3d0
                 An_1(0,j)  = max(An_1(0,j), 0d0)
                 An_1(0,j)  = min(An_1(0,j), 1d0)

             endif

             if (maskC(nx+1,j) .eq. 1) then
             
                 An_1(nx+1,j) = ( 4d0 * An_1(nx,j) - An_1(nx-1,j) ) / 3d0
                 An_1(nx+1,j) = max(An_1(nx+1,j), 0d0)
                 An_1(nx+1,j) = min(An_1(nx+1,j), 1d0)

             endif

         enddo

            
     endif
    
    !-- Looping through the domain --!
    do i = 0, nx+1
        do j = 0, ny+1
            !-- If Ocean --!
            if (maskC(i,j) .eq. 1) then

                !-- Computing the volume fraction --!
                if (mu_phi) then
                    Phi_I(i, j) = max(0d0, Phi_0 - c_phi * inertial(i, j))
                
                else
                    Aout(i, j)= max(0d0, Phi_0 - c_phi * inertial(i, j))

                endif

            endif
        enddo
    enddo

    ! if (peri .ne. 0)   call periodicBC2(Aout)	
        
    return


end subroutine volumefraction_phi


subroutine div(utp, vtp)

    implicit none

    include 'parameter.h'
    include 'CB_mask.h'
    include 'CB_const.h'
    include 'CB_DynVariables.h'

    integer i, j
    double precision, intent(in) :: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
    ! double precision, intent(out):: div(nx,ny)

    do i = 1, nx
        do j = 1, ny

        if (maskC(i,j) .eq. 1) then

            div_I(i,j)=(utp(i+1,j)-utp(i,j) + vtp(i,j+1)-vtp(i,j))/Deltax 

        endif

        enddo
    enddo

end subroutine div

subroutine divergence_muphi()

    use datetime, only: datetime_type
    use ellipse
    use muphi
    
    implicit none

    include 'parameter.h'
    include 'CB_Dyndim.h'
    include 'CB_DynVariables.h'
    include 'CB_DynForcing.h'
    include 'CB_const.h'
    include 'CB_mask.h'
    include 'CB_options.h'

    integer i, j
    double precision tan_dilat_angle, microscopic, pi

    pi = 4d0 * atan( 1d0 )
    microscopic = tan(phi_f_micro*pi/180)

    do i = 0, nx+1
        do j = 0, ny+1

            if (maskC(i, j) .eq. 1) then

                if (Phi_eq) then
                    tan_psi(i,j) = K_div*(A(i,j) - Phi_I(i,j))

                else 
                
                    tan_psi(i, j) = (mu_I(i, j) - microscopic) /(1+microscopic*mu_I(i, j))

                endif

                div_I(i, j) = shearC_I(i, j) * tan_psi(i, j)

            endif
            
        enddo
    enddo

end subroutine divergence_muphi

subroutine SIC_2_PHI

!----------------------------------------------------------!
!--- SUBROUTINE USED TO COMPUTE THE VOLUME FRACTION PHI ---!  
!---              FROM THE EMPIRICAL LAW                ---!    
!----------------------------------------------------------!

    use muphi

    implicit none 

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    !-- Local Variables --!
    integer i, j
    double precision :: diff_A
    

    !-- Looping through the domain --!
    do i = 0, nx+1
        do j = 0, ny+1
            !-- If Ocean --!
            if (maskC(i,j) .eq. 1) then
                diff_A = A(i ,j)-1d0
                Phi_S(i, j) = (1+Phi_min) * A(i,j) + Phi_min
                Phi_G(i, j) = -(1-Phi_max) * tanh((diff_A)/(1d0/2d0)) + (Phi_max)
                Phi_A(i, j) = Phi_S(i, j) * Phi_G(i,j)
            endif
        enddo
    enddo

end subroutine SIC_2_PHI

subroutine non_dimensional_shear

!----------------------------------------------------------!
!--- SUBROUTINE COMPUTING THE NON DIMENSIONAL SHEAR FOR ---!
!---              THE FRICTIONAL PRESSURE               ---!
!----------------------------------------------------------!

    !- Loading modules -!
    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'
    
    !- Local Variables -!
    integer i, j
    double precision eps, highI, lowI, Press

    !-- Looping though the domain --!
    do i = 0, nx+1
        do j = 0, ny+1

            !-- If Ocean at i, j
            if (maskC(i,j) .eq. 1) then

                    ! Capping Pressure 
                    Press = max(Pmax(i, j), 1d-20)

                    !-- Computing I from Pmax --!
                    Ifriction(i, j) = max(SQRT(rhoice * h(i, j)/Press) * d_average*shearC_I(i, j), 1d0)

                    !-- Computing shear rate --!
                    gamma_I(i, j) = c_1*Ifriction(i, j)**c_2

            endif
        enddo
    enddo

    !--- Boundary conditions ---!

    if (Periodic_y .eq. 0) then
        do i = 0, nx+1

            if (maskC(i,0) .eq. 1) then             
                Ifriction(i,1)  = lowI
            endif

            if (maskC(i,ny+1) .eq. 1) then             
                Ifriction(i,ny)  = lowI
            endif

        enddo
    endif            
	  
    if (Periodic_x .eq. 0) then
        do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then             
                Ifriction(1,j)  = lowI
            endif

            if (maskC(nx+1,j) .eq. 1) then             
                Ifriction(nx,j)   = lowI  
            endif           

        enddo
    endif

end subroutine non_dimensional_shear

