!***********************************************************************
!     subroutine Pressure: 
!
!     calculates the ice strenght using:
!       
!       P = Pstar * h * exp(-C (1 - A))
!
!************************************************************************


subroutine Ice_strength ()
   use ellipse
   use muphi

   use, intrinsic :: ieee_arithmetic
   implicit none

   include 'parameter.h'
   include 'CB_DynVariables.h'
   include 'CB_mask.h'
   include 'CB_options.h'
   include 'CB_const.h'


   integer i, j, peri


   peri = Periodic_x + Periodic_y

   if (Rheology .eq. 1) then
      
      !FSTD: changing the domain of calculation
      ! do i = 1, nx
      !    do j = 1, ny
      do i = 0, nx+1
         do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then
               Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
               Pp(i,j) = Pp(i,j) / 2d0
               Pt(i,j) = Tens * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )  
               Pt(i,j) = Pt(i,j) / 2.0d0 
            endif
         enddo
      enddo

   elseif ( Rheology .eq. 2 ) then

!         do i = 1, nx
!            do j = 1, ny
!              if (maskC(i,j) .eq. 1) then
!               Pmax(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ))
!               Pmin(i,j) = Cohe  * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ))
!               p(i,j) = 0.0d0
!              endif
!            enddo
!         enddo

   elseif (Rheology .eq. 4) then 
      do i = 0, nx+1
         do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then
               Pp(i,j) = 2d0 * Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
               Pt(i,j) = 2d0* Tens * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )  
            endif
         enddo
      enddo

   endif

!------- set P = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------

	
   if (Periodic_y .eq. 0) then ! Open or Closed Boundary
      
      do i = 0, nx+1

         if (maskC(i,0) .eq. 1) then  ! Open Boundary           
            Pp(i,1)  = 0d0
            Pt(i,1)  = 0d0
         endif

         if (maskC(i,ny+1) .eq. 1) then             

            !---- Removed the setting P at the open boundary -----!
            !---- Changed the index here ----!
            ! Pp(i,ny+1)  = 0d0
            ! Pt(i,ny+1)  = 0d0    ! Pt(i, 1) replaced by Pt(i, ny)
            Pp(i,ny)  = 0d0
            Pt(i,ny)  = 0d0    ! Pt(i, 1) replaced by Pt(i, ny)
         endif

      enddo
   endif          
	  
   if (Periodic_x .eq. 0) then
      do j = 0, ny+1

         if (maskC(0,j) .eq. 1) then             
            Pp(1,j)  = 0d0
            Pt(1,j)  = 0d0 ! Replaced here
         endif

         if (maskC(nx+1,j) .eq. 1) then             
            Pp(nx,j)  = 0d0  
            Pt(nx,j)  = 0d0 ! Replaced here
         endif           

      enddo
   endif


      ! If periodicBC in x, y or both
      ! if (peri .ne. 0) call periodicBC(Pt, Pp)

   return
end subroutine Ice_strength
      
!************************************************************************
      !
         ! P_min = 1e-3
! double precision eps, P_min, diff_A, weight, A_0, shear, scaled_A, diff_phi
   ! eps = 1d-08
      ! elseif ( Rheology .eq. 4 ) then

      !    if (Pres_sum) then
            
      !       do i = 0, nx+1
      !          do j = 0, ny+1
      !             if (maskC(i,j) .eq. 1) then

      !                diff_A = max(1d-20, 1-A(i, j))
      !                Pmax(i,j) = max(Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) )), 1d-20)

      !                P_col(i, j)  = Pmax(i,j)*tanh((rhoice * h(i, j)* (( d_average * shearC_I(i, j) ) / (diff_A))**2) / Pmax(i,j))

      !                if (A(i, j) <= 0.5) then
      !                   P_fric(i, j) = 0d0
      !                else
      !                   P_fric(i,j) = Pmax(i, j)*dexp(A(i, j)*C*tan_psi(i, j)*gamma_I(i, j))

      !                endif

      !                if ((Pres_f .eqv. .true.) .and. (Pres_c .eqv. .false.)) then
      !                   Pp(i,j) = P_fric(i,j)
                     
      !                elseif ((Pres_c .eqv. .true.) .and. (Pres_f .eqv. .false.)) then

      !                   Pp(i, j) = P_col(i,j)

      !                else
      !                   Pp(i, j) = P_col(i,j) + P_fric(i,j)

      !                endif
                     
      !                Pt(i, j) = 0d0
      !                ! Pp(i, j) = Pmax(i, j)
                     
      !             endif
      !          enddo
      !       enddo

      !    else

      ! !--------- Pressure for mu(I)-Phi(I)          ---------
      ! !          Pmax = Pstar*h*exp(-C(1-A))        ---------
      ! !          Peq = rhoice*h*((d*shear)/(Phi - Phi_0))**2
      ! !          P = min(Peq, Pmax)                 
      !       do i = 0, nx+1
      !          do j = 0, ny+1
      !             if (maskC(i,j) .eq. 1) then

                  
      !                if (A2Phi) then
      !                   diff_A = max(1d-20, 1-Phi_A(i, j))

      !                else 
      !                   diff_A = max(1d-20, 1-A(i, j))
                     
      !                endif

      !                diff_phi = max(1d-20, 1-Phi_I(i, j))

      !                !-- I'm not sure about the Capping on Pmax --!
      !                ! Maximum pressure
      !                Pmax(i,j) = Pstar * h(i,j) * dexp(-C * ( diff_A ))
      !                ! Pdilat(i, j) = Pmax(i, j) * tanh(tan_psi(i, j)*h(i, j))

      !                Peq(i, j) = rhoice * h(i, j) * (( d_average * shearC_I(i, j) ) / (diff_A))**2  
                     
      !                if (P_dilat) then 
      !                   ! Pp(i, j) =  Pstar * h(i,j) * dexp(-C * ( diff_A )-tan_psi(i ,j))
      !                   ! Pp(i, j) =  Pstar * h(i,j) * dexp(-C * ( diff_A )+an_psi(i ,j))
      !                   ! Pp(i, j) = Pmax(i, j)/2*(1-tan_psi(i, j)/abs(tan_psi(i ,j)))
      !                   ! Pp(i, j) = Pmax(i, j)*((h(i,  j)/2)**1/2)
                        
      !                   ! if (h(i,j) < 2) then
      !                   !    Pp(i, j) = (5*(1/2d0)*(9.80)*(rhowater-rhoice)*rhoice/rhowater + &
      !                   !          5/(5-1)*mu_I(i, j)*(rhowater-rhoice)*9.80*(rhoice*(5-1)/rhowater)**2/(2*0.8))*h(i,j)**2/6
      !                   ! else
      !                   !    Pp(i, j) = 0d0
                        
      !                   ! endif
      !                   Pp(i, j) = Pmax(i, j)
      !                ! Pressure from mu phi 

      !                elseif (A2Phi) then

      !                   Pp(i, j) = Pstar*tanh(Peq(i, j)/Pstar)
      !                   ! Pp(i,j) = Peq(i ,j)
                     
      !                else 

      !                   if (regularization .eq. 'capping') then 
                           
      !                      Pp(i, j) = min(Peq(i, j), Pmax(i, j))
                        
      !                   elseif (regularization .eq. 'tanh') then
                           
      !                      if (Pmax(i, j) .lt. 1d-10) then
      !                         Pp(i, j) = 0d0
      !                      else
      !                         Pp(i, j) = Pmax(i, j) * tanh(Peq(i, j) / Pmax(i, j))
      !                         ! Pp(i, j) = Pmax(i, j)

      !                      endif

      !                   endif

      !                   Pt(i, j) = 0d0
      !                   ! Pp(i,j) = Pmax(i, j)

      !                endif
      !             endif
      !          enddo
      !       enddo
      
      !    endif