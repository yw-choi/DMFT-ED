module dmft_lattice
    use constants
    use matsubara_grid
    use numeric_utils
    use dmft_params

    implicit none

    double complex, allocatable :: &
        Hk(:,:,:,:,:,:)  ! Hk(nk,norb,norb,nspin,na,na)
                         ! lattice hamiltonian

contains

    ! @TODO only 2D square lattice hamiltonian is implemented.
    ! This routine needs to be generalized to deal with general lattices
    subroutine setup_lattice_hamiltonian
        integer :: i,j,ik,nkx,iorb,ispin,ia

        double precision :: dk, kx, ky, ckx, cky

        allocate(Hk(nk,norb,norb,nspin,na,na))

        Hk = cmplx(0.0d0, 0.0d0)

        nkx = int(sqrt(float(nk)))
        dk = 2.0D0*pi/float(nkx)

        ik = 0
        do i = 1, nkx
            do j = 1, nkx
                ik = ik + 1
                kx = dk*float(i-1)
                ky = dk*float(j-1)
                ckx = cos(kx)
                cky = cos(ky)
                do ia=1,na
                    do iorb=1,norb
                        do ispin=1,nspin
                            Hk(ik,iorb,iorb,ispin,ia,ia) = -0.5D0*(ckx+cky)
                        enddo
                    enddo
                enddo 
            enddo
        enddo

    end subroutine setup_lattice_hamiltonian

    ! Calculates the lattice Green's function at (iw,ik). 
    ! Only site/orbital/spin diagonal elements are returned.
    ! @TODO needs to be generalized?
    subroutine lattice_green_function(iw, ik, Sigma, Gkout)
        use utils
        integer, intent(in) :: ik, iw
        double complex, intent(in) :: Sigma(nwloc,norb,nspin,na)
        double complex, intent(out) :: Gkout(norb,nspin,na)

        double complex :: invGk(na*norb*nspin,na*norb*nspin),&
            Gk(na*norb*nspin,na*norb*nspin)
        integer :: i,j, ia,iorb, ja,jorb, ispin


        ! invGk = iw + mu - Hk - Sigma
        invGk = cmplx(0.0d0,0.0d0)

        ! @TODO indexing is somewhat arbitrary.
        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    i = (ia-1)*norb*nspin+(ispin-1)*norb+iorb
                    invGk(i,i) = invGk(i,i) + &
                        cmplx(0.0d0,omega(iw)) + mu - Sigma(iw,iorb,ispin,ia)

                    do ja=1,na
                        do jorb=1,norb
                            j = (ja-1)*norb*nspin+(ispin-1)*norb+jorb
                            invGk(j,i) = invGk(j,i) - Hk(ik,iorb,jorb,ispin,ia,ja)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! matrix inversion
        call cinv(invGk,na*norb*nspin,na*norb*nspin,Gk)

        ! Extract diagonal elements to output
        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    i = (ia-1)*norb*nspin+(ispin-1)*norb+iorb
                    Gkout(iorb,ispin,ia) = Gk(i,i)
                enddo
            enddo
        enddo 
    end subroutine lattice_green_function
end module dmft_lattice
