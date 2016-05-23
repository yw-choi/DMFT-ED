module mpi
    implicit none
    include 'mpif.h'

    integer :: &
        nprocs, & ! number of processors
        taskid, & ! rank of the current processor
        mpierr, & ! mpi error code
        comm      ! MPI_Comm_world

    logical :: &
        master    ! taskid.eq.0 

    public
contains
    subroutine init_mpi

        comm = MPI_Comm_world

        call MPI_Init(mpierr)
        call MPI_Comm_size(comm,nprocs,mpierr)
        call MPI_Comm_rank(comm,taskid,mpierr)

        master = taskid.eq.0

    end subroutine init_mpi

    subroutine finalize_mpi
        call MPI_Finalize(mpierr)
    end subroutine finalize_mpi
end module mpi
