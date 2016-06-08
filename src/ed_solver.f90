module ed_solver
    use mpi
    use utils, only: die
    use dmft_params, only: norb, nspin
    use matsubara_grid, only: nwloc
    use ed_projection, only: project_to_impurity_model
    use ed_params, only: ed_read_params, nbath, nsite, nsector, &
                         nev, nstep, sectors, ek_in, vk_in, diag_method, &
                         eigpair_t

    use ed_hamiltonian, only: ek,vk
    use ed_full_diag, only: full_diagonalize

    implicit none

    ! number of eigenvalues that has prob>PROB_THRESHOLD
    ! and corresponding eigenvalues & eigenvectors
    integer :: nev_calc
    type(eigpair_t), allocatable :: eigpairs(:)

contains

    subroutine ed_init
        call ed_read_params

        allocate(ek(nsite,nspin),vk(norb,nbath,nspin))
        ! input levels, hybridization
        ek = ek_in
        vk = vk_in

    end subroutine ed_init

    subroutine ed_solve(G0,Sigma)
        double complex, intent(in) :: G0(nwloc,norb,nspin)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)
        integer :: i,j

        ! Find ek,vk by fitting the cluster quantity(e.g. G0cl) to G0.
        call project_to_impurity_model(G0,ek,vk)

        ! Diagonalize the AIM Hamiltonian characterized by ek,vk and
        ! return the eigpairs.
        select case(diag_method)
            case ("full")
                ! call full_diagonalize(nev_calc,eigpairs)
            case default
                call die("ed_solve", "Diagonalization method is not implemented.")
        end select
    end subroutine ed_solve

end module ed_solver
