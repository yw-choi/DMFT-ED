module ed_operators
    use dmft_params, only: norb
    use ed_params, only: nbath, nsite, kind_basis
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
                basis_out = generate_basis( basis%ne_up+1, basis%ne_down)
            else
                basis_out = generate_basis( basis%ne_up-1, basis%ne_down)
            endif
        else
            if (pm.eq.1) then
                basis_out = generate_basis( basis%ne_up, basis%ne_down+1)
            else
                basis_out = generate_basis( basis%ne_up, basis%ne_down-1)
            endif
        endif

        allocate(vec_out(basis_out%nloc))
        allocate(vec_all(basis%ntot))
        vec_out = 0.0D0

        call mpi_allgatherv(vec,basis%nloc,mpi_double_precision,vec_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,ierr)

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

    integer function get_bitidx(isite,ispin)
        integer :: isite, ispin
        get_bitidx = (ispin-1)*Nsite + isite-1
    end function get_bitidx

    ! number operator eigenvalue. c^+ (isite,ispin) c (isite,ispin)
    integer function num_op(basis_in,isite,ispin)
        integer(kind=kind_basis) :: basis_in
        integer :: isite,ispin
        integer :: bitidx

        bitidx = get_bitidx(isite,ispin)
        if (BTEST(basis_in,bitidx)) then
            num_op = 1
        else
            num_op = 0
        endif
    end function

    subroutine onsite(basis_in,isite,ispin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        double precision, intent(out) :: coeff

        integer :: bitidx,sgn
        bitidx = get_bitidx(isite,ispin)

        call number_op(basis_in,isite,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0D0
            return
        endif
        
        coeff = ek(isite)
    end subroutine onsite

    subroutine hybridization1(basis_in,iorb,ibath,ispin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,ibath,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        double precision, intent(out) :: coeff

        integer :: sgn
        coeff = 1.0D0
        call destruction_op(basis_in,Norb+ibath,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0D0
            return
        endif
        coeff = coeff*sgn
        call creation_op(basis_out,iorb,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0D0
            return
        endif
        coeff = coeff*sgn
        coeff = coeff*vk(iorb,norb+ibath)
    end subroutine hybridization1

    subroutine hybridization2(basis_in,iorb,ibath,ispin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,ibath,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        double precision, intent(out) :: coeff

        integer :: sgn
        coeff = 1.0D0
        call destruction_op(basis_in,iorb,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0D0
            return
        endif
        coeff = coeff*sgn
        call creation_op(basis_out,Norb+ibath,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0D0
            return
        endif
        coeff = coeff*sgn
        coeff = coeff*vk(iorb,norb+ibath)
    end subroutine hybridization2

    subroutine dens_dens(basis_in,iorb,jorb,ispin,jspin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,jorb,ispin,jspin
        integer(kind=kind_basis), intent(out) :: basis_out
        double precision, intent(out) :: coeff

        integer :: n1,n2

        call number_op(basis_in,iorb,ispin,basis_out,n1)
        call number_op(basis_in,jorb,jspin,basis_out,n2)

        if (n1*n2.eq.0) then
            coeff = 0
            return
        endif

        if (iorb.eq.jorb) then
            if (ispin.eq.jspin) then
                coeff = 0.0D0
                basis_out = 0
                return
            endif

            ! intra orbital coupling
            coeff = U
        else if (ispin.eq.jspin) then
            ! inter orbital coupling, same spin
            coeff = U - 3*Jex ! Uprime - J
        else
            ! inter orbital coupling, different spin
            coeff = U - 2*Jex ! Uprime 
        endif

    end subroutine dens_dens

    subroutine spin_flip(basis_in,iorb,jorb,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,jorb
        integer(kind=kind_basis), intent(out) :: basis_out
        double precision, intent(out) :: coeff

        integer :: sgn_tot,sgn

        sgn_tot=1
        coeff = 0.0D0

        if (iorb.eq.jorb) then
            basis_out = 0
            return
        endif

        ! c_j,u
        call destruction_op(basis_in,jorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_j,d
        call creation_op(basis_out,jorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c_i,d
        call destruction_op(basis_out,iorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_i,u
        call creation_op(basis_out,iorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        coeff = -Jex*sgn_tot ! Jprime

        return
    end subroutine spin_flip

    subroutine pair_exchange(basis_in,iorb,jorb,basis_out,)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,jorb
        integer(kind=kind_basis), intent(out) :: basis_out

        integer :: sgn_tot,sgn

        sgn_tot=1

        if (iorb.eq.jorb) then
            basis_out = 0
            return
        endif

        ! c_j,d
        call destruction_op(basis_in,jorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c_j,u
        call destruction_op(basis_out,jorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_i,d
        call creation_op(basis_out,iorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_i,u
        call creation_op(basis_out,iorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        coeff = -Jex*sgn_tot ! Jprime

        return
    end subroutine pair_exchange

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
end module ed_operators
