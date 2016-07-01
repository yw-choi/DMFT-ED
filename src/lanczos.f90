module lanczos

    use numeric_utils, only: mpi_dot_product, mpi_norm

    implicit none

contains

    ! calculates matrix elements of matrix M in the lanczos basis.
    ! M  =  ( a1  b2  0   0   ... 0   ) 
    !       ( b2  a2  b3  0   ... 0   )
    !       ( 0   b3  a3  b4  ... 0   )
    !       ( 0   0   b4  a4  ... 0   )
    !       ( ...             ... bn  )
    !       ( ...             bn  an  )
    !
    ! uses an external routine matmult that handles the matrix-vector product.
    ! refer the interface declared inside the subroutine.
    subroutine lanczos_iteration(matmult, nloc, vec, nstep, a, b)
        integer, intent(in) :: &
            nloc, &       ! dimension of the vector local to the node
            nstep         ! maximum number of iteration steps

        double precision, intent(in) :: &
            vec(nloc)     ! starting vector for lanczos iteration
        double precision, intent(out) :: &
            a(nstep), & ! diagonal matrix element in lanczos basis
            b(nstep)    ! off-diagonal matrix element in lanczos basis

        interface 
            ! external subroutine for matrix multiplication
            ! Y = M*X
            subroutine matmult(n,X,Y)
                integer, intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)
            end subroutine matmult
        end interface

        ! temporary lanczos vectors
        double precision :: v(nloc,2), w(nloc)
        double precision :: norm_v
        
        integer :: j, ierr
        
        ! Lanczos steps
        ! ref: https://en.wikipedia.org/wiki/Lanczos_algorithm#Iteration 

        ! normalize the initial vector
        norm_v = mpi_norm( vec, nloc)
        v(:,2) = vec/norm_v
        v(:,1) = 0.0D0

        a(:) = 0.0D0
        b(:) = 0.0D0

        ! v(:,1) = v_(j-1)
        ! v(:,2) = v_j
        ! w(:)   = w_j
        lanczos_loop: do j=1,nstep-1
            ! w_j = H*v_j
            call matmult( nloc, v(:,2), w(:) )

            ! a_j = dot(w_j,v_j)
            a(j) = mpi_dot_product(w(:), v(:,2), nloc)

            ! w_j = w_j - a_j * v_j - b_j * v_(j-1)
            w(:) = w(:) - a(j)*v(:,2) - b(j)*v(:,1)

            ! b_(j+1) = norm(w_j)
            b(j+1) = mpi_norm(w(:), nloc)

            if (b(j+1).lt.1e-8) then
                exit lanczos_loop
            endif

            ! v_(j-1) = v_j
            v(:,1) = v(:,2)

            ! v_(j+1) = w_j/b(j+1)
            v(:,2) = w(:)/b(j+1)
        enddo lanczos_loop

        ! handles the last step
        call matmult( nloc, v(:,2), w(:) )
        a(nstep) = mpi_dot_product(w(:),v(:,2),nloc)

    end subroutine lanczos_iteration

    ! subroutine lanczos_iteration(matmult,nloc,vec,nstep,a,b)
    !     integer, intent(in) :: &
    !         nloc, &       ! dimension of the vector local to the node
    !         nstep         ! maximum number of iteration steps

    !     double precision, intent(in) :: &
    !         vec(nloc)     ! starting vector for lanczos iteration
    !     double precision, intent(out) :: &
    !         a(nstep), & ! diagonal matrix element in lanczos basis
    !         b(nstep)    ! off-diagonal matrix element in lanczos basis

    !     interface 
    !         ! external subroutine for matrix multiplication
    !         ! Y = M*X
    !         subroutine matmult(n,X,Y)
    !             integer, intent(in) :: n
    !             double precision, intent(in) :: x(n)
    !             double precision, intent(out) :: y(n)
    !         end subroutine matmult
    !     end interface

    !     ! lanczos vectors
    !     double precision, allocatable :: v(:,:), w(:)
    !     double precision :: norm_v
        
    !     integer :: j, ierr
        
    !     allocate(v(nloc,3),w(nloc))

    !     v(1:nloc,2) = 0.0D0
    !     w(1:nloc) = 0.0D0

    !     a(0:nstep) = 0.0D0
    !     b(0:nstep) = 0.0D0

    !     v(:,1) = vec
    !     ! normalize the initial vector
    !     norm_v = mpi_dot_product( v(:,1), v(:,1), nloc)
    !     v(:,1) = v(:,1)/sqrt(norm_v)
    !     call matmult(nloc,v(:,1),w(:))
    !     a(0) = mpi_dot_product(v(:,1),w(:),nloc)
    !     b(0) = 0.0D0

    !     v(:,2) = w(:) - a(0)*v(:,1)
    !     b(1) = mpi_dot_product(v(:,2),v(:,2),nloc)
    !     call matmult(nloc,v(:,2),w(:))
    !     a(1) = mpi_dot_product( w(:), v(:,2), nloc)
    !     a(1) = a(1)/b(1)
    !     v(:,2) = v(:,2)/sqrt(b(1))

    !     lanczos_loop: do j = 2, nstep
    !         v(:,3) = w(:)/sqrt(b(j-1)) - a(j-1)*v(:,2) - sqrt(b(j-1))*v(:,1)

    !         b(j) = mpi_dot_product(v(:,3),v(:,3),nloc)
    !         call matmult(nloc,v(:,3),w(:))

    !         a(j) = mpi_dot_product(v(:,3),w(:),nloc)
    !         a(j) = a(j)/b(j)
    !         v(:,3) = v(:,3)/sqrt(b(j))
    !         if(b(j).lt.1.e-12) then
    !             exit lanczos_loop
    !         endif
    !         v(:,1) = v(:,2)
    !         v(:,2) = v(:,3)
    !     enddo lanczos_loop

    !     deallocate(v,w)
    ! end subroutine lanczos_iteration
end module lanczos
