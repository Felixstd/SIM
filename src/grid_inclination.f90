subroutine grid_inclination

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    integer i, j, x1, x2, y1, y2, row1, row2, col1, col2
    
    double precision :: slope

    slope = 1

    x1 = 0
    y1 = 50

    x2 = 100
    y2 = 0


    do i = 0, nx+1
        row1 = x1 + i 
        col1 = y1 + i * slope
        row2 = x2 + i
        col2 = y2 + i * slope
        if ((row1 < nx) .and. (col1 < nx) .and. (col2 < nx)) then
            do j = col1, col2 + 1
                maskC(i, j) = 1
            enddo
        endif

    enddo

    
end subroutine grid_inclination