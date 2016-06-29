module ed_full_diag
!==============================================================================
! Serial, full diagonalization for each sector using LAPACK.
! This is for testing small problems.
! Currently, only one sector can be diagonalized in a run.
!==============================================================================
    
    use mpi
    use dmft_params, only: norb, beta
    use ed_params, only: nsite, nbath, eigpair_t, sectors, nsector, PROB_THRESHOLD
    use ed_basis, only: generate_basis, basis_t, ed_basis_get
    use ed_hamiltonian, only: generate_hamiltonian, ek,vk
    use numeric_utils, only: boltzmann_factor
    use utils, only: die

    implicit none

    public :: full_diagonalize

    private
contains

    subroutine full_diagonalize(nev_calc,eigpairs)
        integer, intent(out) :: nev_calc
        type(eigpair_t), allocatable, intent(out) :: eigpairs(:)

        ! local variables
        integer :: isector, iev, ne_up, ne_down, nloc, nh
        double precision, allocatable :: H(:,:), ev(:), prob(:)
        double precision :: Z
        type(basis_t) :: basis

        integer :: i,j

        if (nsector.gt.1) then
            call die("full_diagonalize", &
                "full diagonalization for nsector>1 is not implemented.")
        else if (nprocs.gt.1) then
            call die("full_diagonalize", &
                "full diagonalization for nprocs>1 is not implemented.")
        endif

        if (master) then
            write(*,"(x,a,I3)") "Full diagonalization..."
        endif

        ne_up = sectors(1,1)
        ne_down = sectors(1,2)
        nh = sectors(1,3)

        ! basis states for the sector (ne_up,ne_down)
        call generate_basis(ne_up, ne_down, basis)
            
        allocate(H(nh,nh),ev(nh),prob(nh))

        call generate_hamiltonian(basis,H)
        ! open(unit=71,file="ekvk.dump",form="formatted",status="replace")
        ! do i=1,nsite
        !     write(71,"(F20.12)") ek(i,1)
        ! enddo
        ! write(71,*)
        ! do i=1,norb
        !     do j=1,nbath
        !         write(71,"(F20.12)") vk(i,j,1)
        !     enddo
        !     write(71,*)
        ! enddo
        ! close(71)

        ! open(unit=71,file="hamiltonian.dump",form="formatted",status="replace")
        ! do i=1,nh
        !     do j=1,nh
        !         write(71,"(F20.12)", advance="no") H(i,j)
        !     enddo
        !     write(71,*)
        ! enddo
        ! close(71)

        ! open(unit=71,file="basis.dump",form="formatted",status="replace")
        ! do i=1,basis%nloc
        !     write(71,"(I5,5x)",advance="no") i
        !     do j=1,nsite
        !         if (BTEST(ed_basis_get(basis,i),j-1)) then
        !             write(71,"(I1)", advance="no") 1
        !         else
        !             write(71,"(I1)", advance="no") 0
        !         endif
        !     enddo
        !     write(71,*)
        ! enddo
        ! close(71)


        ! note that lapack returned eigenvalues will be in ascending order
        call lapack_diag(basis,H,ev)

        nev_calc = 1
        prob(1) = 1.0d0 ! ground state
        do iev=2,nh
            prob(iev) = boltzmann_factor(beta, ev(iev)-ev(1))
            if (prob(iev).gt.PROB_THRESHOLD) then
                nev_calc = nev_calc + 1
            endif
        enddo
        Z = sum(prob(1:nev_calc))

        allocate(eigpairs(nev_calc))

        do iev=1,nev_calc
            eigpairs(iev)%sector = 1
            eigpairs(iev)%level  = iev
            eigpairs(iev)%idx    = iev
            eigpairs(iev)%val    = ev(iev)
            eigpairs(iev)%prob   = prob(iev)/Z
            allocate(eigpairs(iev)%vec(sectors(1,3)))
            eigpairs(iev)%vec(:) = H(:,iev)
        enddo

        if (master) then
            write(6,*)
            write(6,*) "Obatined eigenvalues."
            write(6,"(a,ES10.3,a,I5)") " Number of eigenvalues (with prob > ", &
                                       PROB_THRESHOLD,") = ", nev_calc
            write(6,*)
            write(6,"(a)") " Eigenvalue          Prob         Sector   Level"
            do iev=1,nev_calc
                write(6,"(1x,F16.10,4x,ES11.5,2I8)") eigpairs(iev)%val,eigpairs(iev)%prob,&
                                          eigpairs(iev)%sector,eigpairs(iev)%level
            enddo
            write(6,*)
        endif

        ! open(unit=71,file="eigenvector.dump",form="formatted",status="replace")
        ! do i=1,nh
        !     do j=1,nev_calc
        !         write(71,"(F20.12)", advance="no") eigpairs(j)%vec(i)
        !     enddo
        !     write(71,*)
        ! enddo
        ! close(71)
        ! stop
    end subroutine full_diagonalize

    subroutine lapack_diag(basis,H,ev)
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: H(basis%ntot,basis%ntot)
        double precision, intent(out) :: ev(basis%ntot)
        ! local variables
        integer :: nh, info
        double precision :: work(3*basis%ntot)

        nh = basis%ntot

        call DSYEV('V','U',basis%ntot,H,basis%ntot,ev,work,3*basis%ntot,info)

        if (info.ne.0) then
            write(*,*) "DSYEV info = ", info
            call die("lapack_diag","failed to diagonalize")
        endif

    end subroutine lapack_diag

end module ed_full_diag
