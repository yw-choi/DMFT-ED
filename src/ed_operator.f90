module ed_operator
    use dmft_params, only: norb
    use ed_params, only: nbath, nsite, kind_basis
    use mpi
    use ed_basis

    implicit none

contains

    ! Apply creation/destruction operator to a state vector.
    ! pm = 1 : create
    ! pm = 2 : destroy
    subroutine apply_c(basis,vec,pm,iorb,ispin,basis_out,vec_out)
        include 'mpif.h'

        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: vec(basis%nloc)
        integer, intent(in) :: pm, iorb, ispin
        type(basis_t), intent(out) :: basis_out
        double precision, allocatable, intent(out) :: vec_out(:)

        double precision, allocatable :: vec_all(:)
        integer(kind=kind_basis) :: basis_i, basis_j
        integer :: i,j,sgn

        if (ispin.eq.1) then
            if (pm.eq.1) then
                call generate_basis( basis%ne_up+1, basis%ne_down, basis_out)
            else
                call generate_basis( basis%ne_up-1, basis%ne_down, basis_out)
            endif
        else
            if (pm.eq.1) then
                call generate_basis( basis%ne_up, basis%ne_down+1, basis_out)
            else
                call generate_basis( basis%ne_up, basis%ne_down-1, basis_out)
            endif
        endif

        allocate(vec_out(basis_out%nloc))
        allocate(vec_all(basis%ntot))
        vec_out = 0.0D0

        call mpi_allgatherv(vec,basis%nloc,mpi_double_precision,vec_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,mpierr)

        do i=1,basis_out%nloc
            basis_i = ed_basis_get(basis_out,i)
            if (pm.eq.1) then
                call destruction_op(basis_i,iorb,ispin,basis_j,sgn)
            else
                call creation_op(basis_i,iorb,ispin,basis_j,sgn)
            endif

            if (sgn.eq.0) then
                cycle
            endif

            j = ed_basis_idx(basis, basis_j)
            vec_out(i) = vec_out(i) + vec_all(j)*sgn
        enddo

        deallocate(vec_all)
    end subroutine apply_c

    integer function sgn(basis_in,i,j)
        integer(kind=kind_basis) :: basis_in
        integer :: i,j
        integer :: k, sgnsum

        do k=i,j
            if (BTEST(basis_in,k)) then
                sgnsum = sgnsum + 1
            endif
        enddo

        if (mod(sgnsum,2)) then
            sgn = -1
        else
            sgn = +1
        endif
        return
    end function

    subroutine destruction_op(basis_in,isite,ispin,basis_out,sgn)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        integer, intent(out) :: sgn

        ! internal variables
        integer :: i, bitidx, n

        bitidx = get_bitidx(isite,ispin)

        ! already occupied ?
        if (BTEST(basis_in,bitidx)) then
            ! how many times the sign is changed.
            n = 0
            do i=0,bitidx-1
                if (BTEST(basis_in,i)) then
                    n = n + 1
                endif
            enddo

            ! set 0 for (isite,ispin)
            basis_out = IBCLR(basis_in,bitidx)
            if (mod(n,2).eq.1) then
                sgn = -1
            else
                sgn = +1
            endif
        else
            ! if the site is not occupied with the same spin, 
            ! the result is the null vector.
            basis_out = 0
            sgn = 0
        end if
    end subroutine destruction_op

    subroutine creation_op(basis_in,isite,ispin,basis_out,sgn)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        integer, intent(out) :: sgn

        ! internal variables
        integer :: i, bitidx, n

        bitidx = get_bitidx(isite,ispin)
        if (.not.BTEST(basis_in,bitidx)) then
            n = 0
            do i=0,bitidx-1
                if (BTEST(basis_in,i)) then
                    n = n + 1
                endif
            enddo
            basis_out = IBSET(basis_in,bitidx)
            if (mod(n,2).eq.1) then
                sgn = -1
            else
                sgn = +1
            endif
        else
            basis_out = 0
            sgn = 0
        end if
    end subroutine creation_op
end module ed_operator
