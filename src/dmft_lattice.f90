module dmft_lattice
    use constants
    use matsubara_grid
    use numeric_utils
    use dmft_params

    implicit none

    double complex, allocatable :: &
        Hk(:,:,:,:)       ! Hk(nspin,na*norb,nk) Lattice Hamiltonian

contains

    ! @TODO only 2D square lattice hamiltonian is implemented.
    ! This routine needs to be generalized to deal with general lattices
    subroutine setup_lattice_hamiltonian
        integer :: i,j,ik,nkx,iorb,ispin

        double precision :: dk, kx, ky

        allocate(Hk(nspin,na*norb,na*norb,nk))

        Hk = cmplx(0.0d0, 0.0d0)

        nkx = int(sqrt(float(nk)))
        dk = 2.0D0*pi/float(nkx)

        ik = 0
        do i = 1, nkx
            do j = 1, nkx
                ik = ik + 1
                kx = dk*float(i-1)
                ky = dk*float(j-1)
                do ispin=1,nspin
                    do iorb = 1,na*norb
                        Hk(ispin,iorb,iorb,ik) = -0.5D0*(cos(kx)+cos(ky))
                    enddo
                enddo
            enddo
        enddo
    end subroutine setup_lattice_hamiltonian

    ! Calculates the lattice Green's function at (ispin,iw,ik). 
    ! Only orbital diagonal elements are returned.
    subroutine lattice_green_function(ispin, iw, ik, Sigma, Gkout)
        use utils
        integer, intent(in) :: ispin, ik, iw
        double complex, intent(in) :: Sigma(nspin,na,norb,nwloc)
        double complex, intent(out) :: Gkout(na,norb)

        double complex :: Gk(na*norb,na*norb),invGk(na*norb,na*norb)
        integer :: i,j, ia,iorb, ja,jorb

        invGk = cmplx(0.0d0, 0.0d0)

        ! First, calculates the inverse of Gk(ispin,iw,ik)
        do i=1,na*norb
            do j=1,na*norb
                if (i.eq.j) then
                    ia = (i-1)/norb+1
                    iorb = mod(i-1,norb)+1
                    ! self-energy is assumed to be diagonal
                    invGk(i,i) = cmplx(0.0d0,omega(iw)) + mu - Sigma(ispin,ia,iorb,iw)
                endif

                invGk(i,j) = invGk(i,j) - Hk(ispin,i,j,ik)
            enddo
        enddo

        ! Invert the matrix to get Gk(ispin,iw,ik)
        call cinv(invGk,na*norb,na*norb,Gk)

        ! Extract diagonal elements to output
        do ia=1,na
            do iorb=1,norb
                i = (ia-1)*norb+iorb
                Gkout(ia,iorb) = Gk(i,i)
            enddo
        enddo 
    end subroutine lattice_green_function
end module dmft_lattice
