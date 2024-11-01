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
    double precision eps, highI, lowI

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

                if (h(i, j) < eps) then
                ! force to highI in water when the height is less than eps
                    inertial(i, j) = lowI
                    ! inertial(i, j) = 0
                else
                    inertial(i, j) = min(SQRT(rhoice * h(i, j)/Pp(i, j)) * d_average*shear_I(i, j), highI)
                endif

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


    ! do i = 1, nx
    !     do j = 1, ny
    do i = 0, nx+1
        do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then

                ! if ( 0 > inertial(i,j)-eps .and. 0 < inertial(i,j)+eps ) then 
                !     mu_I(i, j) = mu_0

                ! elseif (I_0/inertial(i, j) < 1d-02) then
                !     mu_I(i, j) = mu_infty
                
                ! ! else
                ! if (h(i, j) < 1d-6) then
                !     mu_I(i, j) = 0 
                ! else
                !     ! mu_I(i, j) = mu_0 + ( mu_infty - mu_0 ) / ( I_0/inertial(i, j) + 1)
                !     ! mu_I(i, j) = min(mu_0 +  ( mu_infty - mu_0 ) / ( I_0/(1-A(i, j)) + 1 ), mu_infty) 
                !     mu_I(i, j) = max(mu_0 +  ( mu_infty - mu_0 ) / ( I_0/(1-A(i, j)) + 1 ), mu_0) 



                ! mu_I(i, j) = mu_0 +  ( mu_infty - mu_0 ) / ( I_0/(1-A(i, j)) + 1 )


                ! scaled_A = 0.1 + (0.8-0.1)*A(i,j)
                if ((dilatancy .eqv. .true.) .or. (mu_phi .eqv. .false.)) then
                    min_inertial = max(inertial(i, j), 1d-20)

                    if (h(i,j) < 1e-6) then

                        mu_I(i, j) = 0d0

                    else 
                        mu_I(i, j) = mu_0 + ( mu_infty - mu_0 ) / ( I_0/min_inertial + 1) !+ tan_psi(i, j)
                    
                    endif
                    ! mu_b_I(i, j) = 1/(2*mu_I(i,j))
                    mu_b_I(i, j) = mu_b
                
                else
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
    double precision tan_dilat_angle
    double precision K

    K = 2

    do i = 0, nx+1
        do j = 0, ny+1

            if (maskC(i, j) .eq. 1) then
                tan_psi(i, j) = K * (A(i,j) - Phi_I(i, j))

                div_I(i, j) = shear_I(i, j) * tan_dilat_angle

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

    ! write (filename,"('shear_test')")
    ! open(1,file = filename, status = 'unknown')
    ! do j = 0, ny+1
    !     write(1,*) (shear_I(i, j), i = 0, nx+1)
    ! enddo

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