module dmft
    use matsubara_grid
    use impurity_solver, only: solver_init, solve
    use dmft_params
    use dmft_lattice
    use utils
    use mpi
    implicit none

    integer ::  &
        iloop     ! DMFT loop count

    logical :: &
        converged  ! is DMFT loop converged?

    double complex, allocatable :: &
        G_prev(:,:,:,:),   & ! G_prev(nwloc,norb,nspin,na) local Green's function of the previous step
        G(:,:,:,:),        & ! G(nwloc,norb,nspin,na)      local Green's function of the current step
        G0(:,:,:,:),       & ! G0(nwloc,norb,nspin,na)     Weiss field of the current step
        Sigma(:,:,:,:)       ! Sigma(nwloc,norb,nspin,na)  Self-energy of the current step

contains

    subroutine dmft_init
        ! Read general DMFT parameters independent of solver.
        call read_dmft_params

        ! impurity-solver-specific initialization
        if (master) then
            write(*,*) "Initializing impurity solver..."
        endif
        call solver_init

        if (master) then
            write(*,*) "Setting up matsubara grid..."
        endif
        call setup_matsubara_grid


        allocate(G_prev(nwloc,norb,nspin,na)) 
        allocate(G(nwloc,norb,nspin,na))
        allocate(G0(nwloc,norb,nspin,na))
        allocate(Sigma(nwloc,norb,nspin,na))

        G_prev = cmplx(0.0d0,0.0d0)
        G      = cmplx(0.0d0,0.0d0)
        G0     = cmplx(0.0d0,0.0d0)
        Sigma  = cmplx(0.0d0,0.0d0)

        if (master) then
            write(*,*) "Setting up lattice Hamiltonian..."
        endif
        call setup_lattice_hamiltonian

        ! Setting up the initial Weiss field
        call update_local_green_ftn
        call update_weiss_ftn

    end subroutine dmft_init

    subroutine dmft_loop
        integer :: ispin, ia

        if (master) then
            write(6,"(a)") repeat("=",80)
            write(*,*) "Start of DMFT SCF loop."
            write(6,"(a)") repeat("=",80)
        endif

        converged = .false.
        iloop = 0
        do while(.not.converged.and.iloop.le.nloop)
            iloop = iloop + 1
            if (master) then
                write(*,*) "DMFT Loop ", iloop
            endif

            do ia=1,na
                if (master) then
                    write(*,"(1x,a,I1,a,I2,a)") "Solving impurity problem for ia=",ia,")..."
                endif

                call solve(G0(:,:,:,ia), Sigma(:,:,:,ia))
            enddo

            call update_local_green_ftn
            call update_weiss_ftn
            call check_dmft_converged
            call dump_data
        enddo

        if (master) then
            write(6,"(a)") repeat("=",80)
            write(*,*) "End of DMFT SCF loop."
            write(6,"(a)") repeat("=",80)
        endif
    end subroutine dmft_loop

    subroutine dmft_post_processing
        ! @TODO calculates dynamical properties..

    end subroutine dmft_post_processing

    subroutine dmft_finalize
        deallocate(G_prev,G,G0,Sigma)

    end subroutine dmft_finalize

    ! Calculates new local green function from the given self-energy, 
    ! by summing the lattice green function over k.
    subroutine update_local_green_ftn
        integer :: ik, iw, ispin, ia, iorb
        ! Lattice Green's function at (ispin,iw,ik)
        double complex :: Gk(norb,nspin,na)

        if (master) then
            write(*,*) "Updating the local Green's function..."
        endif

        ! Keep the previous Local Green Function
        G_prev = G
        G = cmplx(0.0d0,0.0d0)

        do iw=1,nwloc
            do ik=1,nk
                ! sum Gk to G. take only site/orbital/spin diagonal part
                call lattice_green_function(iw,ik,Sigma,Gk)
                G(iw,:,:,:) = G(iw,:,:,:) + Gk(:,:,:)
            enddo
        enddo

        G = G / nk

    end subroutine update_local_green_ftn

    subroutine update_weiss_ftn

        if (master) then
            write(*,*) "Updating the Weiss field..."
        endif

        G0 = 1/(1/G+Sigma)

    end subroutine update_weiss_ftn

    ! Test DMFT convergence.
    subroutine check_dmft_converged
        integer :: iw, ia, iorb
        double precision :: diff, diffsum

        diff = sum(abs(G_prev(:,:,:,:)-G(:,:,:,:)))
        diff = diff / (na*norb*nspin*nwloc)

        call mpi_allreduce(diff,diffsum,1,mpi_integer,mpi_sum,comm,mpierr)

        if (master) then
            write(*,"(1x,a,I4,a,E12.5)") "scf ",iloop," : diff = ",diffsum
        endif

        if (diffsum < scf_tol) then
            converged = .true.
        endif
    end subroutine check_dmft_converged

    subroutine dump_data
    ! dump current green ftn, weiss ftn, self energy
        use utils
        use io_units
        integer :: iw, iorb, ia, ispin
        character(len=100) :: fn
        if (master) then
            write(fn,"(a5,I3.3,a4)") "dmft.",iloop,".dat"
            open(unit=IO_DEBUG_DUMP,file=fn,status="replace")
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        do iw=1,nwloc
                            write(IO_DEBUG_DUMP,"(7F16.8)") omega(iw), &
                                real(G(iw,iorb,ispin,ia)), aimag(G(iw,iorb,ispin,ia)), &
                                real(G0(iw,iorb,ispin,ia)), aimag(G0(iw,iorb,ispin,ia)), &
                                real(Sigma(iw,iorb,ispin,ia)), aimag(Sigma(iw,iorb,ispin,ia))
                        enddo
                        write(IO_DEBUG_DUMP,*)
                    enddo
                enddo
            enddo
            close(IO_DEBUG_DUMP)
        endif
        
        call mpi_barrier(comm,mpierr)
    end subroutine dump_data
end module dmft
