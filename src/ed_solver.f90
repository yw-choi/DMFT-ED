module ed_solver
    use mpi
    use dmft_params
    use matsubara_grid
    use ed_projection, only: project_to_impurity_model

    implicit none

    integer ::    &
        nbath,    &   ! number of bath sites
        nsite,    &   ! norb+nbath
        nsector,  &   ! number of hamiltonain sector in (Q,Sz) basis
        nev,      &   ! number of the lowest eigenvalues to be computed
        nstep         ! number of steps in calculating continued fraction 

    integer, allocatable :: &
        sectors(:,:)  ! sectors(nsector,2) number of up/down electrons in each sectors

    double precision, allocatable :: &
        ek(:),    &   ! ek(nsite)         impurity/bath onsite energies
        vk(:,:)       ! vk(norb,nbath)  impurity-bath hybridization

contains

    subroutine ed_solve(G0,Sigma)
        double complex, intent(in) :: G0(norb,nwloc)
        double complex, intent(out) :: Sigma(norb,nwloc)
        double complex :: G0_cl(norb,nwloc) ! cluster-projected Weiss field

        integer :: i,j

        ! if (master) then
        !     print *, "Before projection"

        !     print *,ek
        !     print *,vk
        ! endif
        call project_to_impurity_model(G0,ek,vk,norb,nbath,mu,nw)

    end subroutine ed_solve

    subroutine ed_read_options
        use fdf
        use utils
        character(len=200) :: text
        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline

        integer :: i,j

        nbath = fdf_get("DMFT.ED.Nbath",0)
        if (nbath < norb) then
            call die("ed_read_options", "nbath cannot be less than norb.")
        endif
        nsite = norb+nbath

        nsector = fdf_get("DMFT.ED.Nsector",1)
        allocate(sectors(nsector,2))

        if (fdf_block('DMFT.ED.Sectors', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nsector))
                sectors(i,1) = fdf_bintegers(pline,1) ! nup
                sectors(i,2) = fdf_bintegers(pline,2) ! ndown
                i = i + 1
            enddo
        endif

        allocate(ek(nsite),vk(norb,nbath))
        ek(:) = 0.0d0
        vk(:,:) = 0.0d0
        if (fdf_block('DMFT.ED.InitialImpurityLevels', bfdf)) then
            i = 1
            do while( (i .le. norb) .and. (fdf_bline(bfdf, pline)))
                ek(i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        if (fdf_block('DMFT.ED.InitialBathLevels', bfdf)) then
            i = 1
            do while( (i .le. nbath) .and. (fdf_bline(bfdf, pline)))
                ek(Norb+i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        if (fdf_block('DMFT.ED.InitialHybridization', bfdf)) then
            i = 1
            do while( i.le.nbath .and. (fdf_bline(bfdf, pline)))
                do j=1,norb
                    vk(j,i) = fdf_breals(pline,j)
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
    end subroutine ed_read_options

end module ed_solver
