module dmft_green
    use dmft_params
    use dmft_lattice
    use mpi
    use dmft_grid
    use io_units

    implicit none

    double complex, allocatable :: &
        G_prev(:,:,:,:),   & ! G_prev(nwloc,norb,nspin,na)
                             ! local Green's function of the previous step

        G(:,:,:,:),        & ! G(nwloc,norb,nspin,na)
                             ! local Green's function of the current step

        G0(:,:,:,:),       & ! G0(nwloc,norb,nspin,na)
                             ! lattice Weiss field

        Sigma(:,:,:,:)       ! Sigma(nwloc,norb,nspin,na)
                             ! self energy

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
        double complex :: Gk(norb,nspin,na)

        if (master) then
            write(*,*) "Updating the local Green's function..."
        endif

        ! Keep the previous Local Green Function
        G_prev = G
        G = cmplx(0.0d0,0.0d0)

        do iw=1,nwloc
            do ik=1,nk
                ! sum Gk to G. take only site/orbital/spin diagonal part
                call lattice_green_function(iw,ik,Sigma,Gk)
                do ia=1,na
                    do ispin=1,nspin
                        do iorb=1,norb
                            G(iw,iorb,ispin,ia)=G(iw,iorb,ispin,ia)+Gk(iorb,ispin,ia)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        G = G / nk

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

end module dmft_green
