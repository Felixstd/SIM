subroutine conjugate_gradients(b,diagx, &
                                x)

    use nonlocalVP
    use numerical_CG


    implicit none 


    include 'parameter.h'
    include 'CB_DynVariables.h'
    include 'CB_mask.h'
    include 'CB_options.h'
    include 'CB_const.h'

    double precision, dimension(0:nx+2,0:ny+2),intent(in) ::  b
    double precision, dimension(0:nx+2,0:ny+2),intent(in) ::  diagx
    double precision, dimension(0:nx+2,0:ny+2),intent(inout) ::  x



    double precision :: & 
        alpha, &
        beta_cg, &
        norm_res_km1, & 
        norm_res_k, & 
        norm_squared,&
        norm_dq, &
        relative_tol,&
        global_norm, &
        dot_product_cg


    double precision, dimension(0:nx+2,0:ny+2) ::  workspace_x, searchdirs_x,pres_x,res_x
    integer :: k,i,j
    logical :: monitor_cg
    character (len = 5) :: precond_type 

    precond_type = 'diag'
    monitor_cg = .true.


    workspace_x  = 0d0
    searchdirs_x = 0d0
    pres_x       = 0d0
    res_x        = 0d0


    norm_squared = 0d0
    alpha        = 0d0
    beta_cg      = 0d0
    norm_res_km1 = 0d0
    norm_res_k   = 0d0

    call matvec_P(x,workspace_x)

    

    ! print*, workspace_x

    call residual_vec(b,workspace_x,res_x)
    
    ! norm_res_k = dot_product_cg(res_x,res_x)

    call precondition_cg(res_x,diagx,searchdirs_x,precond_type)

    ! print *, 'norm_dq       = ', norm_dq
    ! print *, 'norm_res_k    = ', norm_res_k
    ! print *, 'norm_res_km1  = ', norm_res_km1
    ! print *, 'alpha         = ', alpha
    ! print *, 'beta_cg       = ', beta_cg
    ! print *, 'global_norm   = ', global_norm
    ! print *, 'relative_tol  = ', relative_tol
    ! print *, 'max searchdir = ', maxval(abs(searchdirs_x))
    ! print *, 'max workspace = ', maxval(abs(workspace_x))
    ! print *, 'max res       = ', maxval(abs(res_x))
    ! print *, 'max x         = ', maxval(abs(x))

    norm_res_k = dot_product_cg(res_x, searchdirs_x)

    do k = 0, kmax
        
        global_norm = sqrt(dot_product_cg(res_x,res_x))

        if (k == 0) then 
            relative_tol = global_norm*gamma_cg
        endif 

        if (monitor_cg) then
            print *, 'CG: L2norm is', k,global_norm
        endif

        if (global_norm <= relative_tol) then
            exit
        endif

        call matvec_P(searchdirs_x,workspace_x)

        norm_dq = dot_product_cg(searchdirs_x,workspace_x)

        if ( abs(norm_dq) < 1e-20 ) then
        ! if (monitor_cg) then
        !     ! print *, 'k = ', k
        ! print *, 'norm_dq       = ', norm_dq
        ! print *, 'norm_res_k    = ', norm_res_k
        ! print *, 'norm_res_km1  = ', norm_res_km1
        ! print *, 'alpha         = ', alpha
        ! print *, 'beta_cg       = ', beta_cg
        ! print *, 'global_norm   = ', global_norm
        ! print *, 'relative_tol  = ', relative_tol
        ! print *, 'max searchdir = ', maxval(abs(searchdirs_x))
        ! print *, 'max workspace = ', maxval(abs(workspace_x))
        ! print *, 'max res       = ', maxval(abs(res_x))
        ! print *, 'max x         = ', maxval(abs(x))
        !     print *, 'CG WARNING: breakdown, norm_dq = ', norm_dq
        !     print *, 'CG: restarting from current residual'
        ! endif
        ! restart: reset search direction to current residual
        call precondition_cg(res_x, diagx,searchdirs_x, precond_type)
        cycle
        endif

        alpha = norm_res_k/norm_dq


        ! print *, 'k = ', k
        ! print *, 'norm_dq       = ', norm_dq
        ! print *, 'norm_res_k    = ', norm_res_k
        ! print *, 'norm_res_km1  = ', norm_res_km1
        ! print *, 'alpha         = ', alpha
        ! print *, 'beta_cg       = ', beta_cg
        ! print *, 'global_norm   = ', global_norm
        ! print *, 'relative_tol  = ', relative_tol
        ! print *, 'max searchdir = ', maxval(abs(searchdirs_x))
        ! print *, 'max workspace = ', maxval(abs(workspace_x))
        ! print *, 'max res       = ', maxval(abs(res_x))
        ! print *, 'max x         = ', maxval(abs(x))

        do i = 1,nx
            do j = 1,ny 
                if (maskC(i,j) == 1) then 
                    x(i,j) = x(i,j) + alpha*searchdirs_x(i,j)
                    res_x(i,j) = res_x(i,j) - alpha*workspace_x(i,j)
                endif 
            enddo
        enddo

        norm_res_km1 = norm_res_k

        call precondition_cg(res_x,diagx,pres_x,precond_type)

        norm_res_k = dot_product_cg(res_x,pres_x)

        beta_cg = norm_res_k/norm_res_km1

        do i = 1,nx
            do j = 1,ny 
                if (maskC(i,j) == 1) then 
                    searchdirs_x(i,j) = pres_x(i,j) + beta_cg*searchdirs_x(i,j)
                endif 
            enddo
        enddo


    enddo 



end subroutine conjugate_gradients




subroutine residual_vec(b,    x , & 
                        res_x)

    implicit none 

    include 'parameter.h'
    include 'CB_mask.h'

    double precision, dimension(0:nx+2,0:ny+2),intent(in) ::  b
    double precision, dimension(0:nx+2,0:ny+2),intent(in) ::  x
    double precision, dimension(0:nx+2,0:ny+2),intent(out) ::  res_x
    
    integer :: i,j

    do i = 1,nx 
        do j = 1,ny
            if (maskC(i,j) .eq. 1) then
                res_x(i,j) = b(i,j)-x(i,j)
            endif

        enddo 
    enddo

end subroutine residual_vec


function dot_product_cg(vector1_x,vector2_x)& 
    result(dot)

    implicit none 
    include 'parameter.h'
    include 'CB_mask.h'

    integer :: i,j
    double precision,dimension(0:nx+2,0:ny+2) :: vector1_x,vector2_x
    double precision :: dot, prod(1:nx,1:ny)

    dot = 0d0
    do i = 1,nx
        do j = 1,ny
            if (maskC(i,j) .eq. 1) then 
                prod(i,j) = vector1_x(i,j)*vector2_x(i,j)
                dot = dot+prod(i,j)
                ! print*,prod(i,j)
            endif
        enddo 
    enddo 
    


end function dot_product_cg

subroutine precondition_cg(vx, diagx, wx, &
                         precond_type)

    implicit none 

    include 'parameter.h'
    include 'CB_mask.h'


    double precision,dimension(0:nx+2,0:ny+2),intent(in) :: vx
    double precision,dimension(0:nx+2,0:ny+2),intent(in) :: diagx
    double precision,dimension(0:nx+2,0:ny+2),intent(out):: wx

    integer :: i,j
    character (len=5), intent(in) :: &
        precond_type ! type of preconditioner


    wx = 0d0

    if (precond_type == 'ident') then
        do i = 1, nx
            do j = 1, ny
                if (maskC(i,j) .eq. 1) then
                    wx(i,j) = vx(i,j)
                endif
            enddo
        enddo

    elseif (precond_type == 'diag') then 

        do i = 1, nx
            do j = 1, ny
                if (maskC(i,j) .eq. 1) then
                    if (diagx(i,j) <1d-20) then 
                        wx(i,j) = 0d0
                    else
                        wx(i,j) = vx(i,j)/diagx(i,j)
                    endif
                endif
        enddo
    enddo


    endif

end subroutine precondition_cg