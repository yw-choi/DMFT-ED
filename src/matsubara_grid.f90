module matsubara_grid

    use mpi
    use constants
    use dmft_params

    integer :: &
        nwloc           ! number of matsubara frequencies local to the node
    
    integer, allocatable :: &
        nw_procs(:),  & ! nw_procs(nprocs) number of matsubara frequencies for each processor
        nw_offsets(:)   ! nw_offsets(nprocs) frequency index offsets 

    double precision, allocatable :: &
        omega(:)        ! omega(nwloc) matsubara frequencies local to the node

    public
contains

    subroutine setup_matsubara_grid
        integer :: namw, i
        character(len=200) :: msg

        nwloc = nw/nprocs
        namw = mod(nw,nprocs)
        if(taskid.lt.namw) nwloc=nwloc+1

        allocate(nw_procs(0:nprocs-1),nw_offsets(0:nprocs-1))
        call mpi_allgather(nwloc,1,mpi_integer,nw_procs(0),1,mpi_integer,comm,mpierr)

        nw_offsets(0) = 0
        do i = 1, nprocs-1
            nw_offsets(i) = nw_offsets(i-1) + nw_procs(i-1)
        enddo

        allocate(omega(1:nwloc))
        if(taskid.lt.namw) then
            ishift = taskid*nwloc
            do i = 1, nwloc
                omega(i) = (2.0D0*float(ishift+i-1)+1)*pi/beta
            enddo
        else
            ishift = namw+taskid*nwloc
            do i = 1, nwloc
                omega(i) = (2.0D0*float(ishift+i-1)+1)*pi/beta
            enddo
        endif

        if (master) then
            write(*,'(a)') repeat("=",80)
            write(*,'(4x,a)') "Matsubara Frequency Grid"
            write(*,'(a)') repeat("=",80)
            do i=1,nprocs
                write(*,'(4x,a4,I4,1x,a1,1x,I6)') "node",i,":",nw_procs(i-1)
            enddo
            write(*,'(4x,a,I6)') "Total    : ",nw
            write(*,'(a)') repeat("=",80)
            write(*,*)
        endif

    end subroutine setup_matsubara_grid

end module matsubara_grid
