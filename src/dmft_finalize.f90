subroutine dmft_finalize

    use fdf
    use mpi
    implicit none

    call fdf_shutdown
    call finalize_mpi

end subroutine dmft_finalize
