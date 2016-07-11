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

    integer, public :: &
        nbath,       &   ! number of bath sites
        nbathperorb, &   ! nbath per orbital (when orbital diagonal)
        nsite,       &   ! norb+nbath
        nsector,     &   ! number of hamiltonain sector in (Q,Sz) basis
        nev,         &   ! number of the lowest eigenvalues to be computed
        nstep            ! number of steps in calculating continued fraction 

    integer, allocatable, public :: &
        sectors(:,:)  ! sectors(nsector,3) ne_up/ne_down/nbasis in each sectors

    double precision, allocatable, public :: &
        ek_in(:),    &   ! ek_in(nsite)       initial impurity/bath levels
        vk_in(:,:)       ! vk_in(norb,nbath)  initial impurity-bath hybridization

    character(len=100), public :: diag_method

    logical, public :: &
        print_arpack_stat,  &
        read_h_imp_params,  & ! read initial impurity hamiltonian parameters 
                              ! from file h_imp.params
        em_present, ek_present, vk_present

    type, public :: eigpair_t
        integer :: sector  ! sector index
        integer :: level   ! level index within the sector
        double precision :: val   ! eigenvalue 
        double precision :: prob  ! exp(-beta*val)/Z

        integer :: nloc    ! dimension of the eigenvector local to the processor
        double precision, allocatable :: vec(:) ! eigenvector
    end type eigpair_t

    private
contains

    subroutine ed_read_params
        use fdf
        use utils
        character(len=200) :: text
        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline

        integer :: i,j,ispin

        diag_method = fdf_get("DMFT.ED.Diagonalization", "full") 

        nbath = fdf_get("DMFT.ED.Nbath",0)
        if (mod(nbath,norb)/=0) then
            call die("ed_read_params", "nbath should be multiple of norb.")
        endif
        nbathperorb = nbath/norb
        nsite = norb+nbath
        select case(diag_method)
            case ("full")
                if (nsite>8) then
                    call die("ed_read_params", "full diagonalization for nsite > 8 is not recommended. use arpack.")
                endif
            case default

        end select

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

        read_h_imp_params = fdf_get("DMFT.ED.HamiltonianParamsFromFile", .false.)

        allocate(ek_in(nsite),vk_in(norb,nbath))
        ek_in = 0.0d0
        vk_in = 0.0d0
        em_present = fdf_block('DMFT.ED.InitialImpurityLevels', bfdf)
        if (em_present) then
            i = 1
            do while( (i .le. norb) .and. (fdf_bline(bfdf, pline)))
                ek_in(i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        ek_present = fdf_block('DMFT.ED.InitialBathLevels', bfdf)
        if (fdf_block('DMFT.ED.InitialBathLevels', bfdf)) then
            i = 1
            do while( (i .le. nbath) .and. (fdf_bline(bfdf, pline)))
                ek_in(Norb+i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        vk_present = fdf_block('DMFT.ED.InitialHybridization', bfdf)
        if (vk_present) then
            i = 1
            do while( i.le.nbath .and. (fdf_bline(bfdf, pline)))
                do j=1,norb
                    vk_in(j,i) = fdf_breals(pline,j)
                enddo
                i = i + 1
            enddo
        endif

        nstep = fdf_integer("DMFT.ED.ContinuedFractionStep", 20)
        nev = fdf_integer("DMFT.ED.Nev",20)

        print_arpack_stat = fdf_get("DMFT.ED.ARPACKStat", .false.)

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
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif
    end subroutine ed_read_params
end module ed_params
