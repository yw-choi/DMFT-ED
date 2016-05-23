program dmft_ed
    use impurity_solver, only: solve
    implicit none

    integer :: iloop, nloop
    logical :: converged

    converged=.false.
    nloop = 10

    call dmft_init

    iloop = 0
    dmft_loop: do while(.not.converged.and.iloop<nloop)
        call solve
        call check_dmft_converged(converged)
        call new_local_green_ftn
        call update_weiss_field
        iloop = iloop + 1
    enddo dmft_loop

    call dmft_post_processing

    call dmft_finalize

end program dmft_ed
