! calculates the complete b vector for forcing including air stress + sstilt 
! (R1, R2) with pressure gradient (bu_ind, bv_ind) and terms related to the 
! water drag

subroutine bvect (utp,vtp,rhs)

   use solver_choice
   use basal_param
   use muphi

   implicit none
   
   include 'parameter.h'
   include 'CB_DynVariables.h'
   include 'CB_Dyndim.h'
   include 'CB_const.h'
   include 'CB_mask.h'
   include 'CB_DynForcing.h'
   include 'CB_bathymetry.h'
   include 'CB_options.h'

   integer i, j, peri

   double precision &
                     uwavg(0:nx+2,0:ny+2), & ! uw evaluated at the v-loc
                     vwavg(0:nx+2,0:ny+2)    ! vw evaluated at the u-loc

   double precision speed1p, speed2p, rhs(nvar) 
   double precision uavg, vavg
   double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
   double precision A_at_u, bathy_at_u, CBfactor
   double precision h_at_u, hc, minA, h1,h2, A1,A2, alpha
   double precision v1,v2,va,vb,vc,vd,v_at_u, u_at_v

   minA=0.01d0
   alpha = 1d06
   
   peri = Periodic_x + Periodic_y
   !FSTD removed the following line
   if (peri .ne. 0)   call periodicBC(utp,vtp) ! Periodic boundaries
   
   do j = 1, ny
      do i = 1, nx+1

         if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

            vwavg(i,j) = ( vwatnd(i,j)   + vwatnd(i,j+1) &
                           + vwatnd(i-1,j) + vwatnd(i-1,j+1) ) / 4d0
         endif

      enddo
   enddo

   do j = 1, ny+1
      do i = 1, nx

         if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then
            
            uwavg(i,j) = ( uwatnd(i,j)   + uwatnd(i+1,j) &
                           + uwatnd(i,j-1) + uwatnd(i+1,j-1) ) / 4d0
         endif

      enddo
   enddo


   do j = 1, ny ! u comp
      do i = 1, nx+1

         if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

            vavg = ( vtp(i,j)   + vtp(i,j+1) &
                  + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

            speed1p = sqrt( ( utp(i,j) - uwatnd(i,j) ) ** 2 &
                     + ( vavg - vwavg(i,j)  ) ** 2 )
         
            CdwC1(i,j) = max(Cdw * speed1p, 1d-10)
         
            if (BasalStress) then

               bathy_at_u = min(bathy(i-1,j), bathy(i,j))

               h1=h(i-1,j)
               h2=h(i,j)
               A1=A(i-1,j)
               A2=A(i,j)

               h_at_u = (h1+h2)/2d0 + (h1/2d0)*tanh(alpha*(h1-h2)) + & 
                     (h2/2d0)*tanh(alpha*(h2-h1))

               A_at_u = (A1+A2)/2d0 + (A1/2d0)*tanh(alpha*(A1-A2)) + &
                     (A2/2d0)*tanh(alpha*(A2-A1))

               if ( A_at_u .gt. minA ) then 
                  hc = ( A_at_u * bathy_at_u ) / k1
               else
                  hc = 10000d0
               endif
                                 
               Cbasal1(i,j)=0d0

               if (h_at_u .gt. hc) then

                  va=abs(vtp(i-1,j+1))
                  vb=abs(vtp(i,j+1))
                  vc=abs(vtp(i-1,j))
                  vd=abs(vtp(i,j))
                        
                  v1 = (va+vb)/2d0 - (va/2d0)*tanh(alpha*(va-vb)) - &
                        (vb/2d0)*tanh(alpha*(vb-va))

                  v2 = (vc+vd)/2d0 - (vc/2d0)*tanh(alpha*(vc-vd)) - &
                        (vd/2d0)*tanh(alpha*(vd-vc))
               
                  v_at_u = (v1+v2)/2d0 - (v1/2d0)*tanh(alpha*(v1-v2)) - &
                        (v2/2d0)*tanh(alpha*(v2-v1))

                  speed1p=sqrt( utp(i,j)**2 + v_at_u**2 )

                  Cbfactor=k2/(speed1p+umin)
      !                     Cbasal1(i,j) = Cbfactor * (h_at_u -hc) * dexp(-CC * (1d0 - A_at_u))
                  Cbasal1(i,j)=Cbfactor * (h_at_u - hc) * dexp(-CC * (1d0 - A_at_u))
               else
                  Cbasal1(i,j)=0d0
               endif
                  
            endif
                  
         endif
            
      enddo
   enddo


   do j = 1, ny+1 ! v comp
      do i = 1, nx
            
         if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then

            uavg = ( utp(i,j)   + utp(i+1,j) &
                  + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

            speed2p = sqrt( ( uavg - uwavg(i,j)  ) ** 2 &
                  + ( vtp(i,j) - vwatnd(i,j) ) ** 2 )

            CdwC2(i,j) = max(Cdw * speed2p, 1d-10)
            
            if (BasalStress) then
            
               bathy_at_u = min(bathy(i,j), bathy(i,j-1)) ! in fact at v

               h1=h(i,j-1)
               h2=h(i,j)
               A1=A(i,j-1)
               A2=A(i,j)

               h_at_u = (h1+h2)/2d0 + (h1/2d0)*tanh(alpha*(h1-h2)) + &
                     (h2/2d0)*tanh(alpha*(h2-h1))

               A_at_u = (A1+A2)/2d0 + (A1/2d0)*tanh(alpha*(A1-A2)) + &
                     (A2/2d0)*tanh(alpha*(A2-A1))

               if ( A_at_u .gt. minA ) then
                  hc = ( A_at_u * bathy_at_u ) / k1
               else
                  hc = 10000d0
               endif

               Cbasal2(i,j)=0d0

               if (h_at_u .gt. hc) then

                  va=abs(utp(i,j))
                  vb=abs(utp(i,j-1))
                  vc=abs(utp(i+1,j))
                  vd=abs(utp(i+1,j-1))

                  v1 = (va+vb)/2d0 - (va/2d0)*tanh(alpha*(va-vb)) - &
                        (vb/2d0)*tanh(alpha*(vb-va))

                  v2 = (vc+vd)/2d0 - (vc/2d0)*tanh(alpha*(vc-vd)) - &
                        (vd/2d0)*tanh(alpha*(vd-vc))

                  u_at_v = (v1+v2)/2d0 - (v1/2d0)*tanh(alpha*(v1-v2)) - &
                        (v2/2d0)*tanh(alpha*(v2-v1))

                  speed2p=sqrt( u_at_v**2 + vtp(i,j)**2 )
                  
                  Cbfactor=k2/(speed2p+umin)
                  !                  Cbasal2(i,j) = Cbfactor * (h_at_u -hc) * dexp(-CC * (1d0 - A_at_u))
                     Cbasal2(i,j)=Cbfactor * (h_at_u - hc) * dexp(-CC * (1d0 - A_at_u))
                  else
                     Cbasal2(i,j)=0d0
                  endif
               
            endif
               
         endif

      enddo
   enddo


   do j = 1, ny
      do i = 1, nx+1

         if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

            if (solver .le. 2) then ! Picard or JFNK

               if ((dilatancy .eqv. .true.) .and. (rheology .eq. 4)) then

                     bu(i,j) = bu_ind(i,j) - (P(i,j)+(zetaC(i, j)-etaC(i,j))*shearC_I(i,j)*tan_psi(i,j) - &
                           (P(i-1,j)+(zetaC(i-1, j)-etaC(i-1,j))*shearC_I(i-1,j)*tan_psi(i-1,j) )) / Deltax + & ! P is the replacement pressure 
                              CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                                 vwavg(i,j)  * sintheta_w   )
                  
               else
                  bu(i,j) = bu_ind(i,j) - ( P(i,j) - P(i-1,j) ) / Deltax + & ! P is the replacement pressure 
                     CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                     vwavg(i,j)  * sintheta_w   )

               endif
               

            elseif (solver .eq. 3) then ! EVP solver
 
               bu(i,j) = R1(i,j) + &
                     CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                     vwavg(i,j)  * sintheta_w   )


            endif
               
            !--------------------------------------------
            ! If MEB, put the decohesive terms to rhs :
            !--------------------------------------------
            if ( Rheology .eq. 3) then

               !     d ( sig_xx ) / dx
               bu(i,j) = bu(i,j) + &
                     ( sigxx(i,j)*GammaMEB(i,j) &
                     - sigxx(i-1,j)*GammaMEB(i-1,j) ) &
                     /  Deltax

               !     d ( sig_xy) / dy    B1_2
               bu(i,j) = bu(i,j) + (sigxyB(i,j+1)*GammaMEB_B(i,j+1) - &
                     sigxyB(i,j)*GammaMEB_B(i,j) ) /  Deltax

            endif

         else

            bu(i,j) = 0d0

         endif

      enddo
   enddo

   do j = 1, ny+1
      do i = 1, nx

         if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then

            if (solver .le. 2) then ! Picard or JFNK
                  
               if ((dilatancy .eqv. .true.) .and. (rheology .eq. 4)) then
                  
                  bv(i,j) = bv_ind(i,j) - (P(i,j)+(zetaC(i, j)-etaC(i,j))*shearC_I(i,j)*tan_psi(i,j) - &
                           (P(i,j-1)+(zetaC(i,j-1)-etaC(i,j-1))*shearC_I(i,j-1)*tan_psi(i,j-1)) ) / Deltax + & ! P is the replacement pressure 
                              CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                                 vwavg(i,j)  * sintheta_w   )

               else
                  bv(i,j) = bv_ind(i,j) - ( P(i,j) - P(i,j-1) ) / Deltax + & ! P is the replacement pressure
                     CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                     uwavg(i,j)  * sintheta_w   )
               
               endif


            elseif (solver .eq. 3) then ! EVP solver
               
               bv(i,j) = R2(i,j) + &
                     CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                     uwavg(i,j)  * sintheta_w   )
            

            endif
               
            !--------------------------------------------
            ! If MEB, put the decohesive terms to rhs :
            !--------------------------------------------
            if ( Rheology .eq. 3) then

               !     d ( sig_yy ) / dy
               bv(i,j) = bv(i,j) + &
                     ( sigyy(i,j) * GammaMEB(i,j) &
                     - sigyy(i,j-1) * GammaMEB(i,j-1) ) &
                     /  Deltax 

               !     d ( sig_xy) / dx    B1_2
               bv(i,j) = bv(i,j) + ( sigxyB(i+1,j)*GammaMEB_B(i+1,j) &
                     - sigxyB(i,j)*GammaMEB_B(i,j) ) /  Deltax
               
            endif

         else

            bv(i,j) = 0d0
            
         endif

      enddo
   enddo

!-----------------------------------------------------------------------------
!   Set bu and bv to 0.0 at the 'appropriate' open bc
!-----------------------------------------------------------------------------

   if (Periodic_x .eq. 0) then
   
      do j = 1, ny+1

         bu(1,j)    = 0.0d0   ! Bering Strait
         bu(nx+1,j) = 0.0d0   !

         bv(nx+1,j)    = 0.0d0 

      enddo
   
   endif
    
   if (Periodic_y .eq. 0) then     
   
      do i = 1, nx+1

         bv(i,1)    = 0.0d0   ! North Atlantic
         bv(i,ny+1) = 0.0d0   ! no open bc in current configuration
      
         bu(i,ny+1)    = 0.0d0 

      enddo
   
   endif

      ! ! I think that I should remove this line. 
      ! ! It should not be there since u is already periodic and P should be by default. 
   if (peri .ne. 0)   call periodicBC(bu,bv) ! Periodic boundaries
      ! call variable_post2(bu, bv, 7, 1, 24)

   call transformer (bu,bv,rhs,1)

   return

end subroutine bvect


!-----------------------------------------------------------------------------
! calculates R1 + component of pressure gradient + tendency term un1 part                                                                         
subroutine bvect_ind

    ! Maybe add periodic boundary conditions. on bu_ind. 

   implicit none

   include 'parameter.h'
   include 'CB_DynVariables.h'
   include 'CB_const.h'
   include 'CB_DynForcing.h'
   include 'CB_options.h'

   integer i, j, peri
   double precision hvert

   peri = Periodic_x + Periodic_y

   do j = 1, ny
      do i = 1, nx+1

         hvert = ( h(i,j) + h(i-1,j) ) / 2d0

         bu_ind(i,j) = R1(i,j) - &
                        rhof * hvert * ( vwater(i,j) + vwater(i,j+1) ) / 2d0

         if ( BDF .eq. 0 ) then
            bu_ind(i,j) = bu_ind(i,j) + (rhoice * hvert * un1(i,j)) / Deltat
         
         elseif ( BDF .eq. 1 ) then
            bu_ind(i,j) = bu_ind(i,j) + &
                  (rhoice * hvert * ( 2d0*un1(i,j) - 0.5d0*un2(i,j) ) ) / Deltat
         
         endif

      enddo
   enddo

   do j = 1, ny+1
      do i = 1, nx


         hvert = (h(i,j) + h(i,j-1) ) / 2d0

         bv_ind(i,j) = R2(i,j) + &
                        rhof * hvert * ( uwater(i,j) + uwater(i+1,j) ) / 2d0

         if ( BDF .eq. 0 ) then
            bv_ind(i,j) = bv_ind(i,j) + (rhoice * hvert * vn1(i,j)) / Deltat
         elseif ( BDF .eq. 1 ) then
            bv_ind(i,j) = bv_ind(i,j) + &
                  (rhoice * hvert * ( 2d0*vn1(i,j) - 0.5d0*vn2(i,j) ) ) / Deltat
         endif

      enddo
   enddo

      ! if (peri .ne. 0) call periodicBC(bu_ind, bv_ind)

   return
end subroutine bvect_ind


                  ! bu(i,j) = bu_ind(i,j) - (P(i,j)-(zetaC(i, j)-etaC(i,j))*shearC_I(i,j)*tan_psi(i,j) - &
                           ! (P(i-1,j)-(zetaC(i-1, j)-etaC(i-1,j))*shearC_I(i-1,j)*tan_psi(i-1,j) )) / Deltax + & ! P is the replacement pressure 
                           !    CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                           !       vwavg(i,j)  * sintheta_w   )

                  ! bu(i,j) = bu_ind(i,j) - (P(i,j)-(zetaC(i, j)*shearC_I(i,j)*tan_psi(i,j) - &
                     !          (P(i-1,j)-(zetaC(i-1, j))*shearC_I(i-1,j)*tan_psi(i-1,j) ))) / Deltax + & ! P is the replacement pressure 
                     !             CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                     !                vwavg(i,j)  * sintheta_w   )