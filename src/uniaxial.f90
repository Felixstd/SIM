

subroutine velocityBC_Uniaxial(tstep)

    implicit none

    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_const.h'
    include 'CB_mask.h'
    include 'CB_options.h'

    integer i, j

    integer, intent(in) :: tstep

    double precision :: v_init, a_uniaxial

    a_uniaxial = 5d-04
    v_init = 0d0

    do j =ny, ny+1
        do i = 0, nx
            uice(i, j) = 0d0
            ! vice(i, j) = -1d0*(a_uniaxial*Deltat*tstep + v_init)
            vice(i, j) = -1d0*3d-02
        enddo
    enddo

    j = 1
        do i = 0, nx 
        uice(i, j) = 0d0
        vice(i, j) = 0d0
        enddo
    ! enddo

end subroutine velocityBC_Uniaxial


    