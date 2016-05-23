module dmft
    use matsubara_grid
    use impurity_solver, only: solver_init, solve
    use dmft_params
    implicit none

    integer ::  &
        iloop     ! DMFT loop count

    logical :: &
        converged  ! is DMFT loop converged?

    double complex, allocatable :: &
        G_prev(:,:),  & ! G_prev(norb,nwloc) local Green's function of the previous step
        G(:,:),       & ! G(norb,nwloc)      local Green's function of the current step
        G0(:,:),      & ! G0(norb,nwloc)     Weiss field of the current step
        Sigma(:,:),   & ! Sigma(norb,nwloc)  Self-energy of the current step
        Hk(:,:,:)       ! Hk(norb,norb,nk)   Lattice Hamiltonian

contains

    subroutine dmft_init
        ! Read general DMFT parameters independent of solver.
        call read_dmft_params

        ! impurity-solver-specific initialization
        call solver_init

        call setup_matsubara_grid

        allocate(G_prev(norb,nwloc)) 
        allocate(G(norb,nwloc))
        allocate(G0(norb,nwloc))
        allocate(Sigma(norb,nwloc))

        G_prev(:,:) = cmplx(0.0d0,0.0d0)
        G(:,:) = cmplx(0.0d0,0.0d0)
        G0(:,:) = cmplx(0.0d0,0.0d0)
        Sigma(:,:) = cmplx(0.0d0,0.0d0)

        call setup_lattice_hamiltonian
        call update_local_green_ftn
        call update_weiss_ftn

        call dump_data
    end subroutine dmft_init

    subroutine dmft_loop
        converged=.false.
        iloop = 0
        do while(.not.converged.and.iloop<nloop)
            call solve
            ! call check_dmft_converged(converged)
            ! call new_local_green_ftn
            ! call update_weiss_field
            iloop = iloop + 1
        enddo
    end subroutine dmft_loop

    subroutine dmft_post_processing

    end subroutine dmft_post_processing

    subroutine dmft_finalize

    end subroutine dmft_finalize

    ! @TODO H(k), therefore G is assumed to be diagonal in orbital index
    subroutine update_local_green_ftn
        integer :: ik, iorb, jorb, iw
        double complex :: gk

        ! Keep the previous Local Green Function
        G_prev = G
        
        G(:,:) = cmplx(0.0d0,0.0d0)

        do iw = 1,nwloc
            do ik = 1,nk
                do iorb = 1,norb
                    gk = cmplx(0.0d0,omega(iw))+mu-Hk(iorb,iorb,ik)-Sigma(iorb,iw)
                    G(iorb,iw) = G(iorb,iw) + 1/gk
                enddo
            enddo

            !@TODO why averaging about k instead of summing up?
            G(:,iw) = G(:,iw) / nk
        enddo

    end subroutine update_local_green_ftn

    subroutine update_weiss_ftn

        G0(:,:) = 1/(1/G(:,:)+Sigma(:,:))

    end subroutine update_weiss_ftn

    ! @TODO only square lattice hamiltonian is implemented.
    ! The following routine needs to be refactored out later
    subroutine setup_lattice_hamiltonian
        integer :: i,j,ik,nkx,iorb

        double precision :: dk, kx, ky

        allocate(Hk(norb,norb,nk))

        Hk(:,:,:) = cmplx(0.0d0, 0.0d0)

        nkx = int(sqrt(float(nk)))
        dk = 2.0D0*pi/float(nkx)

        ik = 0
        do i = 1, nkx
            do j = 1, nkx
                ik = ik + 1
                kx = dk*float(i-1)
                ky = dk*float(j-1)
                do iorb = 1,norb
                    Hk(iorb,iorb,ik) = -0.5D0*(cos(kx)+cos(ky))
                enddo
            enddo
        enddo
    end subroutine setup_lattice_hamiltonian

    ! dump current green ftn, weiss ftn, self energy and die
    subroutine dump_data
        use utils
        integer :: iw, iorb
        if (master) then
            open(unit=99,file="dump.dat",status="replace")

            do iorb=1,norb
                do iw=1,nwloc
                    write(99,"(7F16.8)") omega(iw), real(G(iorb,iw)), aimag(G(iorb,iw)), &
                                        real(G0(iorb,iw)), aimag(G0(iorb,iw)), &
                                        real(Sigma(iorb,iw)), aimag(Sigma(iorb,iw))
                enddo
                write(99,*)
            enddo
            call die("dump_data", "data dumped.")
        endif
    end subroutine dump_data
end module dmft
