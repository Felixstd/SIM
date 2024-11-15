PROGRAM PrecisionExample
    ! Use REAL(8) or DOUBLE PRECISION for double precision variables
    REAL(8) :: x, y, sum
    
    x = 1.0D-11  ! Double precision constant
    y = 2.0D-11
    sum = (x + y)*1/2

    PRINT *, "Sum with double precision (REAL(8)):", sum

    ! Using SELECTED_REAL_KIND for even higher precision (if available)
!    INTEGER, PARAMETER :: high_precision = SELECTED_REAL_KIND(15, 307)
 !   REAL(high_precision) :: a, b, high_sum

  !  a = 1.0E-16_high_precision  ! Variable of high precision kind
   ! b = -1.0E-16_high_precision
    !high_sum = a + b

    PRINT *, "Sum with high precision (SELECTED_REAL_KIND):", high_sum
END PROGRAM PrecisionExample
