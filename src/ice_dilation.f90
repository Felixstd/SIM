subroutine  dilatation_angle

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

    if (dilatation .eq. 'constant') then

        do i = 1, nx
                do j = 1, ny
                if (maskC(i,j) .eq. 1) then
                    delta_d(i,j)    = delta_c*deg2rad
                    dilat_coef(i,j) = sin(delta_c*deg2rad)
                endif
            enddo
        enddo

    elseif (dilatation .eq. 'variable') then 

        print*, 'delta(eII) not implemented'
        
        stop
    endif 

end subroutine dilatation_angle 