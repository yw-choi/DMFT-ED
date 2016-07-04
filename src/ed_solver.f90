module ed_solver
    use mpi
    use utils, only: die
    use dmft_params, only: norb, nspin, mu, na
    use matsubara_grid, only: nwloc, omega
    use ed_projection, only: project_to_impurity_model
    use ed_params, only: ed_read_params, nbath, nsite, nsector, &
                         nev, nstep, sectors, ek_in, vk_in, diag_method, &
                         eigpair_t

    use ed_hamiltonian, only: ed_hamiltonian_init, ek,vk
    use ed_green, only: ed_green_init, cluster_green_ftn, G_cl, ap,bp,an,bn

    use ed_diag_full, only: diag_full
    use ed_diag_arpack

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

    subroutine ed_solve(iloop,ia,G0,Sigma)
        integer, intent(in) :: iloop, ia
        double complex, intent(in) :: G0(nwloc,norb,nspin)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)
        integer :: iw,iorb,ispin

        ! Find ek,vk by fitting the cluster quantity(e.g. G0cl) to G0.
        if (iloop>1) call project_to_impurity_model(G0,ek,vk)

        ! Diagonalize the AIM Hamiltonian characterized by ek,vk and
        ! return the eigpairs.
        select case(diag_method)
            case ("full")
                call diag_full(nev_calc,eigpairs)
            case ("arpack")
                call diag_arpack(nev_calc,eigpairs)
            case default
                call die("ed_solve", "Diagonalization method is not implemented.")
        end select

        ! calculates the interacting Green's function using Lanczos method.
        call cluster_green_ftn(ia,nev_calc,eigpairs)

        ! the self-energy 
        call cluster_self_energy(Sigma)
    end subroutine ed_solve

    ! @TODO need to be refactored to more logical structure
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

    subroutine ed_post_processing

        integer :: ia,ispin,iorb,iw,istep, iev

        character(len=100) :: fn

        if (master) then
            open(unit=381,file="lanzos.dat",form="formatted",status="replace")
            write(381,*) na,nspin,norb
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        do iev=1,nev_calc
                            write(381,*) ia,ispin,iorb,iev,nstep
                            write(381,*) (ap(istep,iev,iorb,ispin,ia),istep=1,nstep)
                            write(381,*) (bp(istep,iev,iorb,ispin,ia),istep=1,nstep)
                            write(381,*) (an(istep,iev,iorb,ispin,ia),istep=1,nstep)
                            write(381,*) (bn(istep,iev,iorb,ispin,ia),istep=1,nstep)
                        enddo
                    enddo
                enddo
            enddo
            close(381)


            if (na.ne.1) stop "Not implemented"

            ! @TODO na=1 only for now
            write(*,*) "Calculating spectral functions..."

        endif

        call mpi_barrier(comm,mpierr)

    end subroutine ed_post_processing
end module ed_solver
