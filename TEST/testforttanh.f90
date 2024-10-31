program find_nan
    use, intrinsic :: ieee_arithmetic
    implicit none
    real :: array(6), x
    integer :: i

    ! Initialize the array with some values, including a NaN
    array = [1.0, 2.0, 0.0/1.0, 4.0, 5.0, 0.0]

    print *, 'Array elements: ', array
    x = 0/0
    ! Loop through the array to find NaN values
    do i = 1, size(array)
	
            print *, tanh(array(i))
            print*, x + 1

!        endif
    end do
end program find_nan
