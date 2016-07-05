program main
    use dmft, only: dmft_init, dmft_loop, dmft_post_processing, dmft_finalize
    use fdf, only: fdf_init, fdf_shutdown
    use mpi
    use timer, only: t1_run, t2_run, print_elapsed_time
    implicit none

    call mpi_setup
    call fdf_init('input.fdf', 'fdf.out')

    if (master) then
        write(6,"(a)") repeat("=",80)
        write(6,*) "Multi-orbital DMFT"
        write(6,*) "Number of processors = ", nprocs
        write(6,"(a)") repeat("=",80)
    endif

    t1_run = mpi_wtime(mpierr)
    
    call dmft_init
    call dmft_loop
    call dmft_post_processing
    call dmft_finalize

    t2_run = mpi_wtime(mpierr)

    if (master) then
        call print_elapsed_time("Total Running Time", t1_run, t2_run)
    endif

    call fdf_shutdown
    call mpi_shutdown
end program main
