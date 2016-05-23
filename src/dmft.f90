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
        Sigma(:,:)      ! Sigma(norb,nwloc)  Self-energy of the current step

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

end module dmft
