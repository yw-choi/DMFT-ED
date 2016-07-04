module ed_diag_arpack
!==============================================================================
! Parallel diagonalization using Parallel ARPACK library.
!==============================================================================

    use mpi
    use dmft_params, only: norb, beta
    use ed_params, only: nev, nsite, nbath, eigpair_t, sectors, nsector, &
                         PROB_THRESHOLD
    use ed_basis, only: generate_basis, basis_t, ed_basis_get
    use ed_hamiltonian, only: multiply_H_OTF
    use numeric_utils, only: boltzmann_factor
    use utils, only: die

    implicit none

    public :: diag_arpack

    private

contains

    subroutine diag_arpack(nev_calc, eigpairs)
        integer, intent(in) :: nev_calc
        type(eigpair_t), allocatable, intent(in) :: eigpairs(:)

        type(basis_t) :: basis
        integer :: ne_up, ne_down, nh, nloc, isector
        type(eigpair_t) :: eigpairs_all(nev, nsector)

        do isector=1,nsector
            ne_up = sectors(isector,1)
            ne_down = sectors(isector,2)
            nh = sectors(isector,3)
            
            call generate_basis(ne_up, ne_down, basis)

            call diagonalization(isector, basis, eigpairs_all(:, isector))
        enddo

        ! @TODO sort lowest few eigenvalues 

    end subroutine diag_arpack

    subroutine diagonalization(isector, basis, eig)
        ! include 'debug.h'
        include 'stat.h'
        integer, intent(in) :: isector
        type(basis_t), intent(in) :: basis
        type(eigpair_t), intent(out) :: eig(nev)

        character, parameter :: bmat = 'I'
        character(len=2), parameter :: which = 'SA'

        double precision :: v(basis%nloc, 2*nev), workd(3*basis%nloc), &
                            resid(basis%nloc)

        double precision :: workl(2*nev*(2*nev+8)), d(2*nev,2), sigma, tol

        logical :: select(2*nev)
        integer :: iparam(11), ipntr(11), lworkl, info, ido, nconv, ncv, &
                   maxitr, mode, ishfts, ldv

        double precision :: pdnorm2
        external :: pdnorm2, daxpy

        integer :: i, j, ierr

        ! ndigit = -3
        ! logfil = 6
        ! msaupd = 3

        ldv = basis%nloc
        ncv = 2*nev
        lworkl = ncv*(ncv+8)

        tol = 0.0

        ishfts = 1
        maxitr = 500
        mode   = 1

        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode

        info = 0
        ido = 0

        do
            call pdsaupd( comm, ido, bmat, basis%nloc, which, nev, tol, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if (ido .eq. -1 .or. ido .eq. 1) then
                call multiply_H_OTF(basis, workd(ipntr(1)), workd(ipntr(2)))
            else
                exit
            endif
        enddo

        if ( info .lt. 0 ) then
            print *, ' Error with pdsaupd, info = ', info
            print *, iparam(5)
            call die("diagonalization", "PDSAUPD ERROR")
        else
            call pdseupd(comm, .true., 'All', select, d, v, ldv, sigma, &
                bmat, basis%nloc, which, nev, tol, resid, ncv, v, ldv, &
                iparam, ipntr, workd, workl, lworkl, ierr )
            if ( ierr .ne. 0) then
                print *, ' Error with pdseupd, info = ', ierr
                call die("diagonalization", "pdseupd ERROR")
            else
                nconv =  iparam(5)
                ! do j=1, nconv
                !     call multiply_H_OTF(basis, v(:,j), ax)
                !     call daxpy(basis%nloc, -d(j,1), v(1,j), 1, ax, 1)
                !     d(j,2) = pdnorm2( comm, basis%nloc, ax, 1 )
                ! enddo
                ! call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                !     'Ritz values and direct residuals')
            end if
        endif

        do i=1,nconv
            eig(i)%sector = isector
            eig(i)%level  = i
            eig(i)%idx    = -1 ! not yet determined
            eig(i)%val    = d(i,1)
            allocate(eig(i)%vec(basis%nloc))
            eig(i)%vec    = v(:,i)
            eig(i)%nloc   = basis%nloc
        enddo
    end subroutine diagonalization
end module ed_diag_arpack
