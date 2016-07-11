program main
    use dmft, only: dmft_init, dmft_loop, dmft_post_processing, dmft_finalize
    use fdf, only: fdf_init, fdf_shutdown
    use mpi
    use timer, only: t1_run, t2_run, print_elapsed_time
    implicit none
    
    character(len=100) :: input_file

    call mpi_setup

    if (iargc()==0) then
        input_file = "input.fdf"
    else
        call getarg(1,input_file)
    endif

    if (master) then
        write(6,"(a)") repeat("=",80)
        write(6,*)
        write(6,*) "Multi-orbital DMFT"
        write(6,*) "Number of processors = ", nprocs
        write(6,"(1x,A,A)") "Input file = ", trim(input_file)
        write(6,*)
        write(6,"(a)") repeat("=",80)
    endif

    call fdf_init(input_file, 'fdf.out')

    t1_run = mpi_wtime(mpierr)
    
    call dmft_init
    call dmft_loop
    call dmft_post_processing
    call dmft_finalize

    t2_run = mpi_wtime(mpierr)

    if (master) then
        write(6,"(a)") repeat("=",80)
        write(6,*)
        call print_elapsed_time(" Total Running Time", t1_run, t2_run)
        write(6,*)
        write(6,*) "End of run."
        write(6,*)
        write(6,"(a)") repeat("=",80)
    endif

    call fdf_shutdown
    call mpi_shutdown
end program main
