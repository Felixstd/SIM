subroutine nonlocal_ice_strength()

    use ellipse
    use muphi
    use nonlocalVP

    use, intrinsic :: ieee_arithmetic
    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'


    integer i, j, peri
    double precision :: oneoverl2
    double precision, dimension(0:nx+2,0:ny+2) :: b
    double precision, dimension(0:nx+2,0:ny+2) :: diagx

    oneoverl2 = (1d0/l_coop)**2


   Pp = 2d0*Pp
    
   b = 0d0
   do i = 1,nx
      do j = 1,ny
         if (maskC(i,j) .eq. 1) then 
               b(i,j) = oneoverl2*Pp(i,j)/(max(A(i,j),1d-10)**2d0)
               ! print*, Pp_loc(i,j)
         endif
      enddo
   enddo 
   print*, oneoverl2
   print *, 'max b         = ', maxval(abs(b))
   print *, 'max Pp_loc    = ', maxval(abs(Pp_loc))
   print *, 'min A         = ', minval(A(1:nx,1:ny), mask=maskC(1:nx,1:ny)==1)

   
   call form_diag(diagx)
   call conjugate_gradients(b,diagx, &
                           Pp)

   do i = 0, nx+1
         do j = 0, ny+1
            if (maskC(i,j) .eq. 1) then
               Pp(i,j) = Pp(i,j) / 2d0
            endif
         enddo
   enddo

   ! Pp_loc = 0d0

   if (Periodic_y .eq. 0) then ! Open or Closed Boundary
      
      do i = 0, nx+1

         if (maskC(i,0) .eq. 1) then  ! Open Boundary           
            Pp(i,1)  = 0d0
            Pt(i,1)  = 0d0
         endif

         if (maskC(i,ny+1) .eq. 1) then             

            !---- Removed the setting P at the open boundary -----!
            !---- Changed the index here ----!
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


end subroutine nonlocal_ice_strength



subroutine matvec_P(x, Ax)

   use nonlocalVP


   implicit none 


   include 'parameter.h'
   include 'CB_DynVariables.h'
   include 'CB_mask.h'
   include 'CB_options.h'
   include 'CB_const.h'


   double precision, intent(in) :: x(0:nx+2,0:ny+2)
   double precision, intent(out) :: Ax(0:nx+2,0:ny+2)

   integer i, j
   double precision oneoverl2

   oneoverl2 = (1d0/l_coop)**2

   Ax = 0d0
   
   do j = 1, ny
      do i = 1, nx

         if (maskC(i,j) == 1) then 
               Ax(i,j) = -(x(i+1,j) - 2*x(i,j) + x(i-1,j))/Deltax2 &
                     - (x(i,j+1) - 2*x(i,j) + x(i,j-1))/Deltax2 &
                     +(oneoverl2/(max(A(i,j),1d-10))**2)*x(i,j)
         endif

      enddo
   enddo 


end subroutine matvec_P


subroutine form_diag(diagx)

   use nonlocalVP


   implicit none 


   include 'parameter.h'
   include 'CB_DynVariables.h'
   include 'CB_mask.h'
   include 'CB_options.h'
   include 'CB_const.h'


   double precision, dimension(0:nx+2,0:ny+2):: diagx

   double precision oneoverl2
   integer :: i,j

   oneoverl2 = (1d0/l_coop)**2

      diagx = 0d0
      do i = 1, nx
         do j = 1, ny
            if (maskC(i,j) == 1) then
            diagx(i,j) = oneoverl2 / (max(A(i,j),1d-10))**2

            ! x-direction: each active neighbor contributes +1/Deltax2
            ! to the diagonal (from the -(-2*x(i,j))/Deltax2 term)
            ! if (maskC(i+1,j) == 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2
            ! if (maskC(i-1,j) == 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2

            ! ! y-direction
            ! if (maskC(i,j+1) == 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2
            ! if (maskC(i,j-1) == 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2

            endif


         ! if (maskC(i+1,j) .eq. 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2
         ! if (maskC(i-1,j) .eq. 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2
         ! if (maskC(i,j+1) .eq. 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2
         ! if (maskC(i,j-1) .eq. 1) diagx(i,j) = diagx(i,j) + 1d0/Deltax2

! diagx(i,j) = diagx(i,j) + (oneoverl2/(max(A(i,j),1d-10))**2)
      enddo
   enddo 


end subroutine form_diag