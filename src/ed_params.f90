module ed_params
!===============================================================================
! This module contains the option parameters used only in ed_solver module.
! Note that parameters declared in this module must not be changed after
! the initialization. 
!===============================================================================

    use dmft_params
    use mpi
    use numeric_utils, only: icom

    implicit none

    public :: ed_read_params

    ! Threshold for boltzmann factor 
    double precision, parameter, public :: PROB_THRESHOLD = 0.001D0

    integer, parameter, public :: &
        KIND_BASIS = 4 ! integer kind for basis representation

    integer, public ::    &
        nbath,    &   ! number of bath sites
        nsite,    &   ! norb+nbath
        nsector,  &   ! number of hamiltonain sector in (Q,Sz) basis
        nev,      &   ! number of the lowest eigenvalues to be computed
        nstep         ! number of steps in calculating continued fraction 

    integer, allocatable, public :: &
        sectors(:,:)  ! sectors(nsector,3) ne_up/ne_down/nbasis in each sectors

    double precision, allocatable, public :: &
        ek_in(:,:),    &   ! ek_in(nsite,nspin)       initial impurity/bath levels
        vk_in(:,:,:)       ! vk_in(norb,nbath,nspin)  initial impurity-bath hybridization

    type, public :: eigpair_t
        integer :: sector  ! sector index
        integer :: level   ! level index within the sector
        integer :: idx     ! index of eigval in ascending order
        double precision :: val   ! eigenvalue 
        double precision :: prob  ! exp(-beta*val)/Z

        integer :: nloc    ! dimension of the eigenvector local to the processor
        double precision, allocatable :: vec(:) ! eigenvector
    end type eigpair_t

    character(len=100), public :: diag_method

    private
contains

    subroutine ed_read_params
        use fdf
        use utils
        character(len=200) :: text
        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline

        integer :: i,j

        diag_method = fdf_get("DMFT.ED.Diagonalization", "full")

        nbath = fdf_get("DMFT.ED.Nbath",0)
        if (nbath < norb) then
            call die("ed_read_params", "nbath cannot be less than norb.")
        endif
        nsite = norb+nbath

        nsector = fdf_get("DMFT.ED.Nsector",1)
        allocate(sectors(nsector,3))

        if (fdf_block('DMFT.ED.Sectors', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nsector))
                sectors(i,1) = fdf_bintegers(pline,1) ! nup
                sectors(i,2) = fdf_bintegers(pline,2) ! ndown
                sectors(i,3) = icom(sectors(i,1)+sectors(i,2),sectors(i,1))* &
                               icom(sectors(i,1)+sectors(i,2),sectors(i,2))
                i = i + 1
            enddo
        endif

        allocate(ek_in(nsite,nspin),vk_in(norb,nbath,nspin))
        ek_in = 0.0d0
        vk_in = 0.0d0
        if (fdf_block('DMFT.ED.InitialImpurityLevels', bfdf)) then
            i = 1
            do while( (i .le. norb) .and. (fdf_bline(bfdf, pline)))
                ! @TODO initial imp/bath levels are assumed to 
                ! be independent of spin.
                ek_in(i,1) = fdf_breals(pline,1)
                ek_in(i,:) = ek_in(i,1)
                i = i + 1
            enddo
        endif

        if (fdf_block('DMFT.ED.InitialBathLevels', bfdf)) then
            i = 1
            do while( (i .le. nbath) .and. (fdf_bline(bfdf, pline)))
                ! @TODO initial imp/bath levels are assumed to 
                ! be independent of spin.
                ek_in(Norb+i,1) = fdf_breals(pline,1)
                ek_in(Norb+i,:) = ek_in(Norb+i,1)
                i = i + 1
            enddo
        endif

        if (fdf_block('DMFT.ED.InitialHybridization', bfdf)) then
            i = 1
            do while( i.le.nbath .and. (fdf_bline(bfdf, pline)))
                do j=1,norb
                    ! @TODO initial imp/bath levels are assumed to 
                    ! be independent of spin.
                    vk_in(j,i,1) = fdf_breals(pline,j)
                    vk_in(j,i,:) = vk_in(j,i,1)
                enddo
                i = i + 1
            enddo
        endif

        nstep = fdf_integer("DMFT.ED.ContinuedFractionStep", 20)
        nev = fdf_integer("DMFT.ED.Nev",20)

        if (master) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,'(4x,a)') "Exact Diagoanlization Solver Parameters"
            write(6,'(a)') repeat("=",80)
            write(text,*) "Diagonalization Method"
            write(6,'(3x,a40,2x,a,2x,a)') text, '=', diag_method
            write(text,*) "Number of impurity orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Norb
            write(text,*) "Number of bath orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nbath
            write(text,*) "Number of sites"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nsite
            write(text,*) "Number of eigenvalues to be computed"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', nev
            write(text,*) "Number of steps in continued fraction"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nstep
            !@TODO dump initial conditions (bath levels, hybridization)
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif
    end subroutine ed_read_params


end module ed_params
