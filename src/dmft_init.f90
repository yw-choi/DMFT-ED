subroutine dmft_init
    use mpi
    use fdf
    use utils, only: die
    use dmft, only: read_dmft_options
    use impurity_solver, only: read_solver_options
    implicit none

    call init_mpi
    if (master) then
        write(6,"(a)") repeat("=",80)
        write(6,*) "Multi-orbital DMFT calculation with ED solver"
        write(6,*) "Number of processors = ", nprocs
        write(6,"(a)") repeat("=",80)
    endif
    
    call fdf_init('input.fdf', 'fdf.out')

    ! Read general DMFT parameters independent of solver.
    call read_dmft_options
    ! Read solver-specific input options
    call read_solver_options

end subroutine dmft_init
