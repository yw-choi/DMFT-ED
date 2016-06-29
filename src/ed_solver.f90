module ed_solver
    use mpi
    use utils, only: die
    use dmft_params, only: norb, nspin, mu
    use matsubara_grid, only: nwloc, omega
    use ed_projection, only: project_to_impurity_model
    use ed_params, only: ed_read_params, nbath, nsite, nsector, &
                         nev, nstep, sectors, ek_in, vk_in, diag_method, &
                         eigpair_t

    use ed_hamiltonian, only: ed_hamiltonian_init, ek,vk
    use ed_full_diag, only: full_diagonalize
    use ed_green, only: ed_green_init, cluster_green_ftn, G_cl

    implicit none

    ! number of eigenvalues that has prob>PROB_THRESHOLD
    ! and corresponding eigenvalues & eigenvectors
    integer :: nev_calc
    type(eigpair_t), allocatable :: eigpairs(:)

contains

    subroutine ed_init
        call ed_read_params
        call ed_hamiltonian_init
        call ed_green_init
    end subroutine ed_init

    subroutine ed_solve(iloop,G0,Sigma)
        integer, intent(in) :: iloop
        double complex, intent(in) :: G0(nwloc,norb,nspin)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)
        integer :: iw,iorb,ispin

        character(len=100) :: fn

        write(fn,"(a,I2.2)") "green.dump.",iloop

        ! Find ek,vk by fitting the cluster quantity(e.g. G0cl) to G0.
        ! open(unit=123,file=fn,form="formatted",status="replace")
        ! do iw=1,nwloc
        !     write(123,*) omega(iw), real(G0(iw,1,1)), aimag(G0(iw,1,1))
        ! enddo
        ! close(123)
        call project_to_impurity_model(G0,ek,vk)

        ! Diagonalize the AIM Hamiltonian characterized by ek,vk and
        ! return the eigpairs.
        select case(diag_method)
            case ("full")
                call full_diagonalize(nev_calc,eigpairs)
            case default
                call die("ed_solve", "Diagonalization method is not implemented.")
        end select

        ! calculates the interacting Green's function using Lanczos method.
        call cluster_green_ftn(nev_calc,eigpairs)
        write(fn,"(a,I2.2)") "green.cl.dump.",iloop

        ! Find ek,vk by fitting the cluster quantity(e.g. G0cl) to G0.
        open(unit=123,file=fn,form="formatted",status="replace")
        do iw=1,nwloc
            write(123,*) omega(iw), real(G_cl(iw,1,1)), aimag(G_cl(iw,1,1))
        enddo
        close(123)

        ! the self-energy 
        call cluster_self_energy(Sigma)

    end subroutine ed_solve

    subroutine cluster_self_energy(Sigma)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)
        integer :: iw,iorb,ispin
        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    Sigma(iw,iorb,ispin)=cmplx(0.0d0,omega(iw))+mu &
                        -ek(iorb,ispin)-delta_cl(iw,iorb,ispin)&
                        -1/G_cl(iw,iorb,ispin)
                enddo
            enddo
        enddo
    end subroutine cluster_self_energy

    double complex function delta_cl(iw,iorb,ispin)
        integer :: iw,iorb,ibath,ispin

        delta_cl = cmplx(0.0d0,0.0d0)

        do ibath=1,nbath
            delta_cl = delta_cl + vk(iorb,ibath,ispin)*vk(iorb,ibath,ispin)&
                        /(cmplx(0.0d0,omega(iw))-ek(norb+ibath,ispin))
        enddo

    end function delta_cl
end module ed_solver
