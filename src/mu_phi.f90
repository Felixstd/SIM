! Not sure about the P here
! 
! 
! 
! 
subroutine inertial_number

    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    double precision eps, highI

    !
    ! Need to look for the eps and see if it is the right thing to do. 
    ! Also, for the value of the High I as P -> 0
    eps = 1d-12
    highI = 1d0

    do i = 0, nx+1
        do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then

                inertial(i, j) = min(SQRT(rhoice * h(i, j)/Pp(i, j)) * d_average*shear_I(i, j), highI)

                

                ! if ((0 > shear_I(i,j)-1d-12 .and. 0 < shear_I(i,j)+1d-12) .or. ( 0 > h(i,j)-eps .and. 0 < h(i,j)+eps )) then 
                !     inertial(i, j) = 0d0
                !     ! inertial(i, j) = highI

                ! elseif ((shear_I(i,j) > 1d-15 .and. shear_I(i,j) < -1d-15) .and. (0 > Pp(i,j)-eps .and. 0 < Pp(i,j)+eps)) then
                !     inertial(i, j) = highI

                ! elseif (0 > Pp(i,j)-eps .and. 0 < Pp(i,j)+eps) then
                !     inertial(i, j) = highI

                ! ! elseif ( 0 > shear_I(i,j)-1d-12 .and. 0 < shear_I(i,j)+1d-12 ) then 
                ! !     inertial(i, j) = 0
                ! else
                !     inertial(i, j) = SQRT(rhoice * h(i, j)/Pp(i, j)) * d_average*shear_I(i, j)
                ! endif

            endif
        enddo
    enddo

    if (Periodic_y .eq. 0) then
        do i = 0, nx+1

            if (maskC(i,0) .eq. 1) then             
                inertial(i,1)  = highI
            endif

            if (maskC(i,ny+1) .eq. 1) then             
                inertial(i,ny)  = highI
            endif

        enddo
    endif            
	  
    if (Periodic_x .eq. 0) then
        do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then             
                inertial(1,j)  = highI
            endif

            if (maskC(nx+1,j) .eq. 1) then             
                inertial(nx,j)   = highI  
            endif           

        enddo
    endif

end subroutine inertial_number

subroutine dilatancy
    
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

                ! if ( 0 > h(i,j)-eps .and. 0 < h(i,j)+eps ) then 
                !     Phi_I(i, j) = 0
                ! else
                    Phi_I(i, j) = Phi_0 - c_phi * inertial(i, j)
                
                ! endif
            ! else 
            !     Phi_I(i, j) = 0
            endif


        enddo
    enddo

    ! do j = 1, ny+1 ! periodic in x
    !     Phi_I(0,j)    = Phi_I(nx,j)
    !     Phi_I(nx+1,j) = Phi_I(1,j)
    !     Phi_I(nx+2,j) = Phi_I(2,j)
    ! enddo

    ! do i = 1, nx+1 ! periodic in y
    !     Phi_I(i,0)    = Phi_I(i,ny)
    !     Phi_I(i,ny+1) = Phi_I(i,1)
    !     Phi_I(i,ny+2) = Phi_I(i,2)
    ! enddo


end subroutine dilatancy

subroutine angle_friction_mu 

    use muphi

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    double precision eps

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
                
                ! else
                    mu_I(i, j) = mu_0 + ( mu_infty - mu_0 ) / ( I_0/inertial(i, j) + 1)
                ! endif
            endif
        enddo
    enddo

! 

end subroutine angle_friction_mu


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

    double precision dudx, dvdy, dudy, dvdx, land, lowA
    double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

    ! double precision shear(0:nx+1,0:ny+1)

    land  = -999d0
    shear_I = land
    ! lowA = -888d0
    lowA = 0d0

    do i = 1, nx
        do j = 1, ny

        if (A(i, j) .lt. 0.1) then
        !     shear_I(i, j) = 0d0 
        
        else
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

                  endif

               
            endif
        enddo
    enddo

    j = 0
    do i = 0, nx+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)     = lowA
        endif
    enddo

    j=ny+1
    do i = 0, nx+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)    = lowA
        endif
    enddo

    i=0
    do j=0,ny+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)     = lowA
        endif
    enddo

    i = nx+1
    do j=0,ny+1
        if ( maskC(i,j) .eq. 1 ) then
            shear_I(i,j)     = lowA
        endif
    enddo

    return
end subroutine shear