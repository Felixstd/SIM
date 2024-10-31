subroutine grid_inclination_mask

    use grid_angle

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    
    double precision :: y, x, center_x, center_y
    double precision :: intercept_1, slope
    double precision :: y_end_1, x_end_2, slope_2, intercept_3
!
! Here, intercept_2 > intercept_1
!

    slope = tan(theta)

    center_x = nx / 2.0
    center_y = ny / 2.0

    ! intercept_1 = - d / (2 * cos(theta))
    ! intercept_2 =  d / (2 * cos(theta))

    intercept_1 = intercept_2 - d*SQRT(slope**2 + 1)

    y_end_1 = slope * nx + intercept_1
    x_end_2 = (ny - intercept_2) / slope 
    slope_2 = (y_end_1 - ny)/(nx - x_end_2)
    intercept_3 = y_end_1 - slope_2*nx

    do i = 1, nx+1
        do j = 1, ny+1

            ! x = i - center_x
            ! y = j - center_y
            x = i
            y = j

            if (y < (- x * slope + intercept_2)) then
                maskC(i, j) = 0.0

            elseif (y > ( x * slope_2 + intercept_3)) then
                maskC(i, j) = 0.0
            
            else 
                maskC(i, j) = 1.0

        
            endif
        end do
    end do


    
end subroutine grid_inclination_mask


subroutine grid_inclination_init

    use grid_angle

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

    integer i, j
    
    double precision :: y, x, center_x, center_y
    double precision :: intercept_1, slope

!
! Here, intercept_2 > intercept_1
!
    slope = tan(theta)

    center_x = nx / 2.0
    center_y = ny / 2.0

    ! intercept_2 = d
    intercept_1 = intercept_2 - d*SQRT(slope**2 + 1)

    do i = 1, nx+1
        do j = 1, ny+1
            
            h(i,j)   =  0d0  ! initial ice thick
            A(i,j)   =  0d0  ! initial ice conc

            x = i
            y = j

            if (maskC(i, j) .eq. 1) then
            
                if (y >= (slope * x + intercept_1) .and. y <= ( slope * x + intercept_2)) then
                    A(i, j) = 1d0
                    h(i, j) = 1d0
                endif
            endif
        end do
    end do


end subroutine grid_inclination_init



