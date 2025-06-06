
      subroutine wind_forcing (date, tstep) 
      use datetime, only: datetime_type, datetime_delta_type, operator(+), operator(-)
      use datetime, only: delta_init, datetime_str, seconds, time_set_from_datetime
      use io, only: load_geostrophic_wind, RP!, load_climatological_wind
      use solver_choice
      use muphi
      use UTILS

      implicit none

      include 'parameter.h'
      include 'CB_DynForcing.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_options.h'
      include 'CB_const.h'
      include 'CB_mask.h'

 
      type(datetime_type), intent(in) :: date
      integer, intent(in) :: tstep
      type(datetime_type) :: date1, date2
      type(datetime_delta_type) :: diff, step
      integer :: dseconds, delta
      character(len=2) :: cdelta

      character(LEN=60) file_u, file_v

      integer i, j, imax, jmax 
      integer  ncell, peri
      integer, parameter :: h2sec = 3600
      integer  year, month, day, hour, minute, second, milli

      double precision tauax, tauay, wspeed, rampfactor, Tramp, wspeed_goal
      double precision uairmax, uairmean,alpha, jj
      double precision uair1(0:nx+2,0:ny+2),vair1(0:nx+2,0:ny+2)
      double precision uair2(0:nx+2,0:ny+2),vair2(0:nx+2,0:ny+2)
      double precision uuair(0:nx+2,0:ny+2),uvair(0:nx+2,0:ny+2)
 
      logical :: verbose = .false.

      peri = Periodic_x + Periodic_y
     
      delta =  int(deltax)/1000

      year = date%year
      month = date%month
      day = date%day
      hour = date%hour
      minute = date%minute
      second = date%second
      milli = date%milli

!------------------------------------------------------------------------
!     load wind data for a given date and time
!------------------------------------------------------------------------

      if ( Wind .eq. '6hours' ) then

         ! Winds are available at 00:00, 06:00, 12:00 and 18:00 hours. 
         step = delta_init(hours=6)
         dseconds = int(modulo(seconds(time_set_from_datetime(date)), seconds(step)))
         
         select case(dseconds)
            
       case(0)
          write(*,*) "Loading geostrophic winds: ", datetime_str(date)
          call load_geostrophic_wind(date, delta, uair, vair)
          
       case default
          
          ! Load winds before and after date. 
          diff = delta_init(seconds=dseconds)
          date1 = date-diff
          date2 = date1+step
          write(*,*) "Loading geostrophic winds: ", datetime_str(date1)
          write(*,*) "Loading geostrophic winds: ", datetime_str(date2)
          call load_geostrophic_wind(date1, delta, uair1, vair1)
          call load_geostrophic_wind(date2, delta, uair2, vair2)
          
          ! Interpolate
          alpha = dseconds/seconds(step)
          uair = alpha *(uair2 - uair1) + uair1
          vair = alpha *(vair2 - vair1) + vair1
          
       end select

       if (verbose) then
          write(*,*) "   dseconds: ", dseconds
          write(*,*) "   alpha:    ", alpha  
       end if
!------------------------------------------------------------------------------

      elseif ( Wind .eq. '60yrs_clim' ) then 
         WRITE(cdelta, '(I2)') delta
         file_u = "forcing/wind/" // cdelta // "/climatological_uwnd_" // cdelta // ".txt"
         file_v = "forcing/wind/" // cdelta // "/climatological_vwnd_" // cdelta // ".txt"         

         open ( unit = 20, file = file_u, status = 'old' )
         
         print *, 'opening/reading ', file_u

         read(20,*)                                               &
              ( ( uair(i,j), i = 1, nx+1 ), j = 1, ny+1 )
         
         close(20)


            
         open ( unit = 21, file = file_v, status = 'old' )
          
         print *, 'opening/reading ', file_v
         
         read(21,*)                                               &
              ( ( vair(i,j), i = 1, nx+1 ), j = 1, ny+1 )
         
         close(21)
         
      elseif ( Wind .eq. 'specified' ) then
         
         wspeed = 10.0d0
         Tramp=6d0*3600d0

         if (RampupWind) then
            wspeed = 10.0d0
            Tramp = 6d0*3600d0

            rampfactor=1d0-exp(-1d0*tstep*Deltat/Tramp)


            if ((rheology .eq. 4) .or. (rheology .eq. 1)) then

            ! linear ramp up 
      

               wspeed = 10d0
               ! Tramp = 60*60*10
               Tramp = 60*60*2
               rampfactor = tanh(tstep*Deltat/Tramp)
               ! rampfactor =  1

            else
               wspeed = 10.0d0
               Tramp = 6d0*3600d0

               rampfactor=1d0-exp(-1d0*tstep*Deltat/Tramp)

            endif

         elseif (RampupForcing) then
            wspeed = 14.0d0
            Tramp = 0.3d0/3.6d3  ! Rate of forcing increase, in N/m2 per hour
            rampfactor = (Tramp * (milli*1d-3 + 1d0*second + &
                                  6.0d1*minute + 3.6d3*hour + &
                                  8.64d4*(day-1d0))  / Cda   )**0.5d0 / wspeed
            rampfactor = min(rampfactor, 1d0)                        
         else
            wspeed = 10.0d0
            ! wspeed = 0.0d0
            rampfactor=1d0
         endif

         print *, 'Specified uair(m/s) =', tstep, rampfactor*wspeed

         do i = 1, nx+1
            jj = 0d0
            do j = 1, ny+1

               if (inclined) then
                  vair(i, j) = -wspeed*rampfactor*cos(theta)
                  uair(i, j) = -wspeed*rampfactor*sin(theta)

               elseif (uniaxial) then

                  if (ny == 1000) then
                     if (j .gt. 500) then
                        vair(i, j) = -wspeed*rampfactor
                        uair(i, j) = 0d0
                     endif
                  else 
                     vair(i, j) = -wspeed*rampfactor
                     uair(i, j) = 0d0
                  endif
               
               elseif (shear_test) then
                  if (j .gt. 100) then 
                     ! uair(i, j) = wspeed*rampfactor
                     ! vair(i, j) = 0d0
                     uair(i, j) = wspeed*rampfactor*cos(theta*PI/180d0)
                     vair(i, j) = -wspeed**rampfactor*sin(theta*PI/180d0)
                  ! !    
                  endif
                  ! uair(i, j) = wspeed*rampfactor*(1-abs(1-2*jj/dble(ny)))

                  ! ! uair(i, j) = wspeed*rampfactor*jj/dble(ny)
                  ! vair(i, j) = 0d0
                  ! uair(i, j) = wspeed*rampfactor*cos(theta)
                  ! vair(i, j) = -wspeed**rampfactor*sin(theta)
               
               else 
                  uair(i,j) = rampfactor*wspeed
                  vair(i,j) = 0d0

               endif

               jj = jj +1d0
               
!     call random_number(rdnumb)
!               uair(i,j) = wm * (rdnumb - 0.5d0)
!               call random_number(rdnumb)
!     vair(i,j) = wm * (rdnumb - 0.5d0)
               
            enddo
         enddo
         
      endif
      
!------------------------------------------------------------------------
!     Calculate wind speed and wind stress at the grid nodes
!     ?speediw is never calculated?
!------------------------------------------------------------------------


      do i = 1, nx+1
         do j = 1, ny+1
            
            speeda(i,j) = sqrt ( uair(i,j) ** 2 + vair(i,j) ** 2 )
            
            uuair(i,j)  = speeda(i,j) * uair(i,j)
            uvair(i,j)  = speeda(i,j) * vair(i,j)
               
         enddo
      enddo

      
      do i = 1, nx+1
         do j = 1, ny+1
            
            
!------------------------------------------------------------------------
!     Calculate the wind curl B-grid. 
!     Used in the reduced gravity ocean model.
!------------------------------------------------------------------------

             
!            curlua(i,j) = (vair(i+1,j) - vair(i-1,j)) / (2d0 * Deltax)
!     +                  - (uair(i,j+1) - uair(i,j-1)) / (2d0 * Deltax)


!------------------------------------------------------------------------
!     Calculate the forcing (tau_a and sea surface tilt) on the C and
!     B-grid, (R1, R2 and R1n and R2n, respectively).
!------------------------------------------------------------------------

            tauax = (                                             &
                      Cda  * ( uuair(i,j)   * costheta_a          &
                             - uvair(i,j)   * sintheta_a )        &
                    + Cda  * ( uuair(i,j+1) * costheta_a          &
                             - uvair(i,j+1) * sintheta_a )        &
                    ) / 2d0

            tauay = (                                             &
                      Cda  * ( uvair(i,j)   * costheta_a          &
                                    + uuair(i,j)   * sintheta_a ) &
                    + Cda  * ( uvair(i+1,j) * costheta_a          &
                                    + uuair(i+1,j) * sintheta_a ) &
                    ) / 2d0

            if (solver .le. 2) then ! Picard or JFNK  

               R1(i,j) = tauax ! 2nd term calc in bvect_ind because of IMEX
               R2(i,j) = tauay ! 2nd term calc in bvect_ind because of IMEX

            elseif (solver .eq. 3) then ! EVP solver
            
               R1(i,j) = tauax &
                       - rhof * ( h(i,j) + h(i-1,j) ) / 2d0   &
                       * ( vwater(i,j)   + vwater(i,j+1) ) / 2d0

               R2(i,j) = tauay &
                       + rhof * ( h(i,j) + h(i,j-1) ) / 2d0   &
                       * ( uwater(i,j)   + uwater(i+1,j) ) / 2d0

            endif

!------------------------------------------------------------------------
!     The sea surface tilt term (rho f h (k x u_w^g) is not included in 
!     the forcing onthe B-grid, since the freedrift solution is solved
!     for (u_i - u_w^g)
!------------------------------------------------------------------------


            R1n(i,j) = Cda   * ( uuair(i,j)   * costheta_a &
                               - uvair(i,j)   * sintheta_a )

            R2n(i,j) = Cda   * ( uvair(i,j)   * costheta_a &
                               + uuair(i,j)   * sintheta_a )



            if ( j .eq. ny+1 ) R1(i,j) = 0d0
            if ( i .eq. nx+1 ) R2(i,j) = 0d0

         enddo
      enddo

      uairmax = 0d0
      uairmean = 0d0
      imax = 0
      jmax = 0
      ncell = 0

      do j=1,ny+1
         do i=1,nx+1

            if ( maskC(i,j) .eq. 1 ) then

               uairmean = uairmean + speeda(i,j)
               ncell = ncell + 1

               if (speeda(i,j) .gt. uairmax) then
                  uairmax = speeda(i,j)
                  imax = i
                  jmax = j
               endif

            endif

         enddo
      enddo

      print *, 'uair max (km/h)=  ', uairmax*3.6d0, '  at  ', imax, jmax
      print *, 'uair mean(km/h)=  ', uairmean*3.6d0/ncell

      if (uairmax*3.6d0 .gt. 300d0) then

         print *, 'WARNING WARNING UAIRMAX WARNING WARNING'
         stop

      endif

      ! if (peri .ne. 0) then
      !    call periodicBC(R1, R2)
      !    call periodicBC(R1n, R2n)
      !    call periodicBC(uair, vair)
      !    call periodicBC(uuair, uvair)
      ! endif

      return
      end
      
