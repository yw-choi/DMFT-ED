module dmft_green
    use dmft_params
    use dmft_lattice, only: dos, dosw, hk
    use mpi
    use dmft_grid
    use io_units
    use numeric_utils, only: cinv

    implicit none

    public :: &
        dmft_green_init, &
        update_local_green_ftn, &
        update_weiss_ftn, &
        dump_green

    double complex, allocatable, public :: &
        G_prev(:,:,:,:),   & ! G_prev(nwloc,norb,nspin,na)
                             ! local Green's function of the previous step

        G(:,:,:,:),        & ! G(nwloc,norb,nspin,na)
                             ! local Green's function of the current step

        G0(:,:,:,:),       & ! G0(nwloc,norb,nspin,na)
                             ! lattice Weiss field

        Sigma(:,:,:,:)       ! Sigma(nwloc,norb,nspin,na)
                             ! self energy

    private 
contains

    subroutine dmft_green_init
        allocate(G_prev(nwloc,norb,nspin,na))
        allocate(G(nwloc,norb,nspin,na))
        allocate(G0(nwloc,norb,nspin,na))
        allocate(Sigma(nwloc,norb,nspin,na))

        G_prev = cmplx(0.0d0,0.0d0)
        G      = cmplx(0.0d0,0.0d0)
        G0     = cmplx(0.0d0,0.0d0)
        Sigma  = cmplx(0.0d0,0.0d0)

        ! Setting up the initial Weiss field
        call update_local_green_ftn
        call update_weiss_ftn

    end subroutine dmft_green_init

    ! Calculates new local green function from the given self-energy,
    ! by summing the lattice green function over k.
    subroutine update_local_green_ftn
        integer :: ik, iw, ispin, ia, iorb
        ! Lattice Green's function at (ispin,iw,ik)
        double complex :: gk(norb,nspin,na), gf

        if (master) then
            write(*,*) "Updating the local Green's function..."
        endif

        ! Keep the previous Local Green Function
        G_prev = G

        ! calculate new local Green's function from the lattice Green's function
        G = cmplx(0.0d0,0.0d0)
        if (tbham==0) then
            ! general tight-binding hamiltonian.
            ! obtain the lattice green's function through matrix inversion.
            do iw=1,nwloc
                do ik=1,nk
                    call gkdiag(iw, ik, sigma, gk)
                    G(iw,:,:,:) = G(iw,:,:,:) + gk(:,:,:)
                enddo
            enddo
            G = G / nk
        else if (tbham==1) then
            ! degenerate bands on 2d square lattice
            ! tight-binding hamiltonian is diagonal. no need of matrix inv.
            do iw=1,nwloc
                do ia=1,na
                    do ispin=1,nspin
                        do iorb=1,norb
                            do ik=1,nk
                                gf = cmplx(0.d0,omega(iw))+mu &
                                    -hk(ik,iorb,iorb,ispin,ia,ia) &
                                    -sigma(iw,iorb,ispin,ia)
                                gf = 1.d0/gf
                                G(iw,iorb,ispin,ia) = G(iw,iorb,ispin,ia)+gf
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            G = G / nk
        else if (tbham==2) then
            ! Bethe lattice.
            ! rather than sum over k, integrate over E using dos
            do iw=1,nwloc
                do ispin=1,nspin
                    do iorb=1,norb
                        call green_bethe(iw,ispin,iorb,gf) 
                        G(iw,iorb,ispin,1) = gf
                    enddo
                enddo
            enddo
        else
            ! other cases should be prevented in the initialization step
            stop "This should not happen"
        endif

    end subroutine update_local_green_ftn

    subroutine update_weiss_ftn

        if (master) then
            write(*,*) "Updating the Weiss field..."
        endif

        G0 = 1.d0/(1.d0/G+Sigma)

    end subroutine update_weiss_ftn

    subroutine dump_green
        ! dump current green ftn, weiss ftn, self energy
        integer :: iw, iorb, ia, ispin
        character(len=100) :: fn
        double precision :: zq
        double complex :: G_all(nw), G0_all(nw), Sigma_all(nw)

        if (master) then
            open(unit=IO_GR_DATA, file="green.dat", status="replace")
            write(IO_GR_DATA, *) na,nspin,norb,nw
        endif

        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    call mpi_allgatherv(G(:,iorb,ispin,ia), &
                                        nwloc, &
                                        mpi_double_precision, &
                                        G_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_precision, &
                                        comm, &
                                        mpierr)
                    call mpi_allgatherv(G0(:,iorb,ispin,ia), &
                                        nwloc, &
                                        mpi_double_precision, &
                                        G0_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_precision, &
                                        comm, &
                                        mpierr)
                    call mpi_allgatherv(Sigma(:,iorb,ispin,ia), &
                                        nwloc, &
                                        mpi_double_precision, &
                                        Sigma_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_precision, &
                                        comm, &
                                        mpierr)

                    if (master) then
                        write(IO_GR_DATA,*) ia, ispin, iorb
                        do iw=1,nw
                            write(IO_GR_DATA,"(7F16.8)") omega(iw), &
                                real(G_all(iw)), &
                                aimag(G_all(iw)), &
                                real(G0_all(iw)), &
                                aimag(G0_all(iw)), &
                                real(Sigma_all(iw)), &
                                aimag(Sigma_all(iw))
                        enddo
                    endif

                    call mpi_barrier(comm,mpierr)
                enddo
            enddo
        enddo

        if (master) then
            close(IO_GR_DATA)
        endif

        call mpi_barrier(comm,mpierr)
    end subroutine dump_green

    ! diagonal element of G_lat(iw,ik)
    subroutine gkdiag(iw,ik,sigma,gk)
        integer, intent(in) :: iw, ik
        double complex, intent(in) :: sigma(nwloc, norb, nspin, na)
        double complex, intent(out) :: gk(norb, nspin, na)

        double complex :: invGk(na*norb*nspin, na*norb*nspin),&
            Gkmat(na*norb*nspin,na*norb*nspin)
        integer :: i,j, ia,iorb, ja,jorb, ispin

        ! invGk = iw + mu - Hk - Sigma
        invGk = cmplx(0.0d0,0.0d0)

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
        call cinv(invGk, na*norb*nspin, na*norb*nspin, Gkmat)

        ! Extract diagonal elements to output
        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    i = (ia-1)*norb*nspin+(ispin-1)*norb+iorb
                    gk(iorb,ispin,ia) = Gkmat(i,i)
                enddo
            enddo
        enddo
    end subroutine gkdiag
    
    subroutine green_bethe(iw,ispin,iorb,gf)
        integer, intent(in) :: iw, ispin, iorb
        double complex, intent(out) :: gf
        integer :: ie
        double precision :: e

        gf = cmplx(0.d0,0.d0)
        ! trapezoidal rule integration
        ! note that at E=-1,1, DOS(E) is exactly 0 so that
        ! we do not care about the boundary term
        do ie=2,nwmodel-1
            e = dosw(ie)

            gf = gf + dos(ie,iorb,ispin,1)/(cmplx(0.d0,omega(iw))+mu-e &
                                    -sigma(iw,iorb,ispin,1))
        enddo

        gf = gf*2.d0/(nwmodel-1)

    end subroutine green_bethe

end module dmft_green
