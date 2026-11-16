subroutine micro_friction()

    use MC_teardrop

    implicit none 

    include 'parameter.h'
    include 'CB_options.h'
    include 'CB_mask.h'
    include 'CB_DynVariables.h'

    integer i, j 
    double precision deg2rad, pi 

    pi = 4 * datan(1d0)

    deg2rad = pi/180d0

    if (friction .eq. 'constant') then

        do i = 1, nx
                do j = 1, ny
                if (maskC(i,j) .eq. 1) then
                    phi_micro(i,j) = phi_prime*deg2rad
                endif
            enddo
        enddo

    elseif (friction .eq. 'muofI') then 

        print*, 'Mu(I): not implemented yet'
        
        stop
    endif 


end subroutine micro_friction 


subroutine macro_friction(phi, delta, mu)

    use MC_teardrop

    implicit none 

    include 'parameter.h'
    include 'CB_options.h'
    include 'CB_mask.h'
    include 'CB_DynVariables.h'

    double precision, intent(in) :: phi(0:nx+2,0:ny+2), delta(0:nx+2,0:ny+2)
    double precision, intent(out) :: mu(0:nx+2,0:ny+2)

    double precision :: phi_macro 
    integer :: i, j 


    do i = 1, nx
        do j = 1, ny
            if (maskC(i,j) .eq. 1) then
                    phi_macro = phi(i,j) + delta(i,j)
                    mu(i,j) = sin(phi_macro)
            endif
        enddo
    enddo


end subroutine macro_friction

