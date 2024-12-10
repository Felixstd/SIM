! Not sure about the P here

subroutine inertial_number

    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    double precision eps, highI, lowI, Press

    !
    ! Need to look for the eps and see if it is the right thing to do. 
    ! Also, for the value of the High I as P -> 0
    eps = 1d-06
    ! highI = 1d0
    highI = 1d0
    lowI  = 0d0

    do i = 0, nx+1
        do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then

                    Press = max(Pp(i, j), 1d-20)
                    inertial(i, j) = min(SQRT(rhoice * h(i, j)/Pp(i, j)) * d_average*shear_I(i, j), highI)

            endif
        enddo
    enddo

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

subroutine non_dimensional_shear

!
! This is to compute the non-dimensional shear gamma^*.
!

    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    double precision eps, highI, lowI, Press

    do i = 0, nx+1
        do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then

                    Press = max(Pmax(i, j), 1d-20)
                    Ifriction(i, j) = max(SQRT(rhoice * h(i, j)/Press) * d_average*shear_I(i, j), 1d0)
                    gamma_I(i, j) = c_1*Ifriction(i, j)**c_2

            endif
        enddo
    enddo

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

subroutine volumefraction_phi
    
    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    double precision eps

    eps = 1d-08

    do i = 0, nx+1
        do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then

                Phi_I(i, j) = max(0d0, Phi_0 - c_phi * inertial(i, j))

            endif


        enddo
    enddo



end subroutine volumefraction_phi

subroutine angle_friction_mu 

    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    double precision eps, scaled_A, diff_A, min_inertial

    eps = 1d-12

    do i = 0, nx+1
        do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then

                ! if ((dilatancy .eqv. .true.) .and. (mu_phi .eqv. .false.)) then
                !     min_inertial = max(inertial(i, j), 1d-20)

                    if (h(i,j) < 1e-2) then
                    ! if (h(i,j) < 1e-6) then

                !     !     !UNTIL 23
                            mu_I(i, j) = 0d0
                            

                    else 
                !     ! mu_I(i, j) = mu_0 + ( mu_infty - mu_0 ) / ( I_0/min_inertial + 1)
                !     mu_I(i, j) = max(mu_0 + ( mu_infty - mu_0 ) / ( I_0/min_inertial + 1) + tan_psi(i, j), 0d0)
                    
                
                        diff_A = max(1d-20, 1-A(i, j))
                        mu_I(i, j) = mu_0 +  ( mu_infty - mu_0 ) / ( (I_0*c_phi)/diff_A + 1 )

                    endif

            endif
        enddo
    enddo

end subroutine angle_friction_mu

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

                div_I(i, j) = shear_I(i, j) * tan_psi(i, j)

            endif
            
        enddo
    enddo

end subroutine divergence_muphi


subroutine shear(utp, vtp)

    use datetime, only: datetime_type
    use ellipse
    
    implicit none

    include 'parameter.h'
    include 'CB_Dyndim.h'
    include 'CB_DynVariables.h'
    include 'CB_DynForcing.h'
    include 'CB_const.h'
    include 'CB_mask.h'
    include 'CB_options.h'

    integer i, j
    character filename*64

    double precision dudx, dvdy, dudy, dvdx, land, lowA
    double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

    ! double precision shear(0:nx+1,0:ny+1)

    land  = -999d0
    shear_I = land
    ! lowA = -888d0
    lowA = 0d0

    do i = 0, nx+1
        do j = 0, ny+1
        ! if (A(i, j) .lt. 0.1) then
        ! if (h(i, j) .lt. 1d-06) then
            ! shear_I(i, j) = 0d0
        
        ! else
            dudx       = 0d0
            dvdy       = 0d0
            dudy       = 0d0
            dvdx       = 0d0
               
            if ( maskC(i,j) .eq. 1 ) then
                  
                    dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                    dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax
                  
                if     ( maskC(i+1,j) + maskC(i-1,j) .eq. 2 ) then
                     
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
                    
                    dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                            ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                            ( 4d0 * Deltax )
                     
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


                shear_I(i,j) = sqrt(( dudx - dvdy )**2d0 &  
                       + ( dudy + dvdx )**2d0 )
                ! div_I(i, j) = (dudx + dvdy)
            
            else
                shear_I(i,j) = land
                ! div_I(i, j) = land

            endif


               
            ! endif
        enddo
    enddo

    j = 0
    do i = 0, nx+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)     = lowA
            ! div_I(i, j) = lowA
        endif
    enddo

    j=ny+1
    do i = 0, nx+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)    = lowA
            ! div_I(i, j) = lowA
        endif
    enddo

    i=0
    do j=0,ny+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)     = lowA
            ! div_I(i, j) = lowA
        endif
    enddo

    i = nx+1
    do j=0,ny+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)     = lowA
            ! div_I(i, j) = lowA
        endif
    enddo


    return
end subroutine shear


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