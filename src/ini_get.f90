!***********************************************************************
!     Subroutine ini_get: set the initial conditions
!***********************************************************************
subroutine ini_get (restart, expno_r, restart_date)
    use datetime, only: datetime_type
    use muphi
    implicit none

    include 'parameter.h'
    include 'CB_options.h'
    include 'CB_DynVariables.h'
    include 'CB_ThermoVariables.h'
    include 'CB_ThermoForcing.h'
    include 'CB_mask.h'
    include 'CB_DynForcing.h'
    include 'CB_const.h'
    include 'CB_buoys.h'
    include 'CB_Dyndim.h'

    character(LEN=32) filename  ! restart file name

    integer, intent(in) :: restart, expno_r
    integer :: i, j, k, year, month, day, hour, minute
    double precision :: Length
    TYPE(datetime_type), intent(in) :: restart_date

    year = restart_date%year
    month = restart_date%month
    day = restart_date%day
    hour = restart_date%hour
    minute = restart_date%minute

!------------------------------------------------------------------------
!     Set initial conditions for h,A,u,v,Ta,Ti,Tl
!------------------------------------------------------------------------

    if ( restart .eq. 0 ) then

       do i = 0, nx+1
          do j = 0, ny+1               
!
!     already zeros on edges

             h(i,j)   =  1d0 * maskC(i,j) ! initial ice thick
             A(i,j)   =  1d0 * maskC(i,j) ! initial ice conc
            !  A(i,j)   =  0.90 * maskC(i,j) ! initial ice conc

!     h and A set to zero on open boundaries

             if ((i.eq.0 .or. i.eq.nx+1) .and. Periodic_x .eq. 0) then ! Changed for Periodic_x
                h(i,j) = 0d0
                A(i,j) = 0d0   
             endif
             if ((j.eq.0 .or. j.eq.ny+1) .and. Periodic_y .eq. 0) then ! Changed for Periodic_y
                h(i,j) = 0d0
                A(i,j) = 0d0  
             endif


!     Uniaxial loading experiment: set bands of open water at the top and sides
            if (Water_Col) then

               if (step_water) then

                  if ((nx == 100) .and. (ny == 250)) then
                     if (i .lt. 21 .or. i .gt. 80) h(i,j) = 0d0
                     if (i .lt. 21 .or. i .gt. 80) A(i,j) = 0d0
                     if (j .gt. 250) h(i,j) = 0d0
                     if (j .gt. 250) A(i,j) = 0d0 


                  elseif ((nx == 200) .and. (ny == 500)) then
                     if (i .lt. 21 .or. i .gt. 180 ) h(i,j) = 0d0
                     if (i .lt. 21 .or. i .gt. 180) A(i,j) = 0d0
                     if (j .gt. 500) h(i,j) = 0d0
                     if (j .gt. 500) A(i,j) = 0d0 

                  
                  elseif ((nx == 200) .and. (ny == 1000)) then
                     if (i .lt. 21 .or. i .gt. 180 ) h(i,j) = 0d0
                     if (i .lt. 21 .or. i .gt. 180) A(i,j) = 0d0

                     !  if (i .lt. 11 .or. i .gt. 100 ) h(i,j) = 0d0
                     !  if (i .lt. 11 .or. i .gt. 190) A(i,j) = 0d0
                     !  if (j .gt. 1000) h(i,j) = 0d0
                     !  if (j .gt. 500) A(i,j) = 0d0 

                  elseif ((nx == 400) .and. (ny == 1000)) then
                     if (i .lt. 101 .or. i .gt. 300) h(i,j) = 0d0
                     if (i .lt. 101 .or. i .gt. 300) A(i,j) = 0d0
                     if (j .gt. 500) h(i,j) = 10d0
                     if (j .gt. 500) A(i,j) = 1d0 

                  elseif ((nx == 200) .and. (ny == 600)) then
                     if (i .lt. 21 .or. i .gt. 180 ) h(i,j) = 0d0
                     if (i .lt. 21 .or. i .gt. 180) A(i,j) = 0d0
                     if (j .gt. 300) h(i,j) = 10d0
                     if (j .gt. 300) A(i,j) = 1d0 

                  elseif ((nx == 102) .and. (ny == 402)) then
                     if (i .lt. 21 .or. i .gt. 81) h(i,j) = 0d0
                     if (i .lt. 21 .or. i .gt. 81 ) A(i,j) = 0d0
                     if (j .lt. 1) h(i,j) = 0d0
                     if (j .lt. 1) A(i,j) = 0d0 

                  endif
               
               else 

                  if ((nx == 200) .and. (ny == 500)) then

                     Length = nx * Deltax / 3
                     A(i, j) = dexp(-(abs(i - nx/2)*Deltax/Length)**6)
                     h(i, j) = dexp(-(abs(i - nx/2)*Deltax/Length)**6)
                     ! if (i .lt. 21 .or. i .gt. 180 ) h(i,j) = 0d0
                     ! if (i .lt. 21 .or. i .gt. 180) A(i,j) = 0d0
                     if (j .gt. 500) h(i,j) = 0d0
                     if (j .gt. 500) A(i,j) = 0d0 

                  endif
               
               endif


            
            endif

             Pp(i,j) = 0d0
             Pt(i,j) = 0d0 
             P(i,j)  = 0d0 

             etaC(i,j)= 0d0
             zetaC(i,j) = 0d0
             etaB(i,j)  = 0d0

             GammaMEB(i,j) = 1d0
             GammaMEB_B(i,j) = 1d0
             dam(i,j) = 1d0
             damB(i,j) = 1d0
             dfactor(i,j) = 1d0
             dfactorB(i,j) = 1d0

             Ta(i,j)  =  273.15d0              ! air temp

             To(i,j)  =  Tof * maskC(i,j)      ! ocean temp
             Ti(i,j)  =  273.15d0 * maskC(i,j)    ! jfl ice surface temp
             Tl(i,j)  =  1d0                   ! land surface temp
               
          enddo
       enddo
       
      !  print*, inclined
       if (inclined) then
         ! print*, 'here'
         call grid_inclination_init
      endif
         
       uice = 0d0 
       vice = 0d0
       un1  = 0d0  ! u-velocity pts 
       vn1  = 0d0  ! v-velocity pts 
       hn1  = h
       An1  = A

    endif

!------------------------------------------------------------------------
!     Load restart files for h,A,u,v,Ta,Ti,Tl initial conditions
!------------------------------------------------------------------------

    if ( restart .eq. 1 ) then               ! load restart files

!------------------------------------------------------------------------
!     Open file and assign unit to it 
!------------------------------------------------------------------------
         
       if ( Dynamic ) then

          write (filename,'("output/h",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (16, file = filename, status = 'old')
            
          write (filename,'("output/A",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (17, file = filename, status = 'old')
            
          write (filename,'("output/u",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (18, file = filename, status = 'old')
            
          write (filename,'("output/v",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (19, file = filename, status = 'unknown')

         ! if ( Rheology .eq. 4) then
         !    write (filename,'("output/p",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
         !          year, month, day, hour, minute, expno_r
         !    open (20, file = filename, status = 'unknown')
            
         !    do j = 1, ny+1
         !       read (20,*) ( P(i,j), i = 1, nx )
         !    enddo

         ! endif 
          
          do j = 1, ny
             read (18,*) ( uice(i,j), i = 1, nx+1 )
          enddo
            
          do j = 1, ny+1
             read (19,*) ( vice(i,j), i = 1, nx )
          enddo

         
            
          do j = 0, ny+1
               
             read (16,*) ( h(i,j),    i = 0, nx+1)
             read (17,*) ( A(i,j),    i = 0, nx+1)
               
          enddo
            
          do k = 16, 19
             close(k)
          enddo
            
          un1 = uice
          vn1 = vice
          hn1 = h
          An1 = A

       endif
         
       Cbasal1 = 0d0
       Cbasal2 = 0d0
         
       if ( Thermodyn ) then
            
          write (filename,'("output/Ta",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (21, file = filename, status = 'unknown')
            
          write (filename,'("output/To",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r            
          open (22, file = filename, status = 'unknown')

          write (filename,'("output/Ti",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (24, file = filename, status = 'unknown')
          
          write (filename,'("output/Qoa",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (73, file = filename, status = 'unknown')
            
          write (filename,'("output/Qsh_io",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (74, file = filename, status = 'unknown')
            
          write (filename,'("output/Pvap",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
               year, month, day, hour, minute, expno_r
          open (75, file = filename, status = 'unknown')
            
          do j = 1, ny
             read (24,*) ( Ti(i,j),   i = 1, nx )
          enddo
            
            
          do j = 0, ny+1
             read (21,*) ( Ta(i,j), i = 0, nx+1 )
             read (22,*) ( To(i,j), i = 0, nx+1 )
             read (73,*) ( Qoa(i,j), i = 0, nx+1 )
             read (74,*) ( Qsh_io(i,j), i = 0, nx+1 )
             read (75,*) ( Pvap(i,j), i = 0, nx+1 )
          enddo
            
          do k = 21, 25
             close(k)
          enddo
            
          do k = 73, 75
             close(k)
          enddo
          
       endif
         
    endif
      
      
  end subroutine ini_get


subroutine initial_conditions_muPhi 

   use ellipse 

   implicit none

   include 'parameter.h'
   include 'CB_options.h'
   include 'CB_DynVariables.h'
   include 'CB_mask.h'
   include 'CB_DynForcing.h'
   include 'CB_const.h'
   include 'CB_buoys.h'

   integer i, j 

   do i = 1, nx
      do j = 1, ny
         if (maskC(i,j) .eq. 1) then
            Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
            ! Phi_I(i,j) = A(i,j)
         ! else 
            ! Pp(i,j) = 0
         endif
      enddo
   enddo

end subroutine initial_conditions_muPhi


subroutine initial_conditions_uniaxial 

   use muphi 

   implicit none

   include 'parameter.h'
   include 'CB_options.h'
   include 'CB_DynVariables.h'
   include 'CB_mask.h'
   include 'CB_DynForcing.h'
   include 'CB_const.h'
   include 'CB_buoys.h'

   integer i, j, c

   ! print *, 'here'

   do i = 0, nx+1
      do j = 0, ny+1

            if (h(i, j) < 1e-6) then

               !BEFORE:
                  ! mu_I(i, j) = mu_infty
               
               !EXPERIMENT 23
               ! mu_I(i, j) = 0d0

               mu_I(i, j) = mu_0
               Phi_I(i, j) = 0d0
               ! Phi_I(i, j) = 0.5
               ! inertial(i, j) = 1d0
               ! inertial(i, j) = 1d0
               inertial(i, j) = 0d0

            ! if (h(i, j) .eq. 0d0) then 
               
            
            else
               ! mu_I(i, j) = 0d0
                mu_I(i, j) = mu_0
               ! Phi_I(i, j) = 1d0
               Phi_I(i, j) = 1
               inertial(i, j) = 1d-16
               
               endif

      enddo
   enddo

end subroutine initial_conditions_uniaxial
