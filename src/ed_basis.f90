module ed_basis

    use mpi
    use dmft_params, only: norb
    use ed_params, only: nbath, KIND_BASIS, nsite
    use numeric_utils, only: icom
    
    implicit none

    public :: generate_basis
    public :: ed_basis_get
    public :: ed_basis_idx
    public :: get_bitidx

    type, public :: basis_t
        integer :: nloc

        integer :: ntot
        integer :: nup
        integer :: ndown

        integer :: ne_up
        integer :: ne_down

        ! For up,down basis, 4-bit integer is sufficient,
        ! because the maximum number of sites will not be more than 31.
        integer, allocatable :: up(:)
        integer, allocatable :: down(:)

        integer, allocatable :: idx_up(:)
        integer, allocatable :: idx_down(:)

        integer, allocatable :: nlocals(:)
        integer, allocatable :: offsets(:)
    end type basis_t

    private
contains

    subroutine generate_basis( ne_up, ne_down, basis )
        integer, intent(in) :: ne_up, ne_down
        type(basis_t), intent(out) :: basis

        ! local variables
        integer :: ispin, i, j, nam
        integer :: minrange, maxrange, counts, nbit
        integer :: nud(2)

        basis%ne_up = ne_up
        basis%ne_down = ne_down

        basis%nup   = icom(nsite,ne_up)
        basis%ndown = icom(nsite,ne_down)
        basis%ntot  = basis%nup * basis%ndown

        basis%nloc  = basis%ntot/nprocs
        nam = mod(basis%ntot,nprocs)
        if (taskid.lt.nam) basis%nloc = basis%nloc + 1

        allocate(basis%nlocals(0:nprocs-1),basis%offsets(0:nprocs-1))
        call mpi_allgather(basis%nloc,1,mpi_integer,basis%nlocals(0),1,mpi_integer,comm,mpierr)

        basis%offsets(0) = 0 
        do i = 1, nprocs-1
            basis%offsets(i) = basis%offsets(i-1) + basis%nlocals(i-1)
        enddo

        allocate(basis%up(basis%nup), basis%down(basis%ndown))

        nud(1) = ne_up
        nud(2) = ne_down

        do ispin=1,2
            ! ref : arXiv:1307.7542, Appendix A
            minrange = 0
            maxrange = 0

            do i=1,nud(ispin)
                minrange = minrange + 2**(i-1)
                maxrange = maxrange + 2**(nsite-i)
            enddo
            
            if (ispin.eq.1) then
                allocate(basis%idx_up(minrange:maxrange))
            else
                allocate(basis%idx_down(minrange:maxrange))
            endif
            
            counts = 0

            do i=minrange,maxrange
                nbit = 0
                do j=0,nsite-1
                    if (BTEST(i,j)) then
                        nbit = nbit + 1
                    endif
                enddo

                if (nbit.eq.nud(ispin)) then
                    counts = counts + 1
                    if (ispin.eq.1) then
                        basis%up(counts) = i
                        basis%idx_up(i) = counts
                    else
                        basis%down(counts) = i
                        basis%idx_down(i) = counts
                    endif
                endif
            enddo
        enddo
    end subroutine generate_basis

    ! ref : arXiv:1307.7542 eq (6)
    integer(kind=kind_basis) function ed_basis_get(basis,idx_loc) 
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: idx_loc

        ! local variables
        integer :: iup, idown, idx

        ! local idx to global idx
        idx = idx_loc + basis%offsets(taskid)

        iup = mod(idx-1,basis%nup)+1
        idown = (idx-1)/basis%nup+1

        ed_basis_get = basis%up(iup)+2**(nsite)*basis%down(idown)
    end function ed_basis_get

    ! ref : arXiv:1307.7542 eq (8)
    integer function ed_basis_idx(basis, basis_i)
        type(basis_t), intent(in) :: basis
        integer(kind=kind_basis) :: basis_i

        ! local variables
        integer :: basis_i_up, basis_i_down
        
        basis_i_up = mod(basis_i,2**(nsite))
        basis_i_down = basis_i/(2**(nsite))

        ed_basis_idx = (basis%idx_down(basis_i_down)-1)*basis%nup + &
                             basis%idx_up(basis_i_up)
    end function ed_basis_idx

    integer function get_bitidx(isite,ispin)
        integer :: isite, ispin
        get_bitidx = (ispin-1)*Nsite + isite-1
    end function get_bitidx
end module ed_basis
