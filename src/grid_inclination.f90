subroutine grid_inclination

    use grid_angle


    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j
    
    double precision :: y, x, center_x, center_y
    double precision :: intercept_1, intercept_2, slope

!
! Here, intercept_2 > intercept_1
!
    slope = tan(theta)

    center_x = nx / 2.0
    center_y = ny / 2.0

    intercept_1 = - d / (2 * cos(theta))
    intercept_2 =  d / (2 * cos(theta))

    


    do i = 1, nx+1
        do j = 1, ny+1

            x = i - center_x
            y = j - center_y

            if (y < -intercept_2 * slope * x / intercept_1 + intercept_2) then
                maskC(i, j) = 0.0
            
            elseif (y >= slope * x + intercept_1 .and. y <= slope * x + intercept_2) then
                maskC(i, j) = 1.0
        
            endif
        end do
    end do


    
end subroutine grid_inclination