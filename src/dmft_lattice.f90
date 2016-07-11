module dmft_lattice
    use mpi
    use constants
    use dmft_grid
    use dmft_params
    use utils, only: die
    use io_units 

    implicit none

    public :: &
        dmft_lattice_init

    double complex, allocatable, public :: &
        Hk(:,:,:,:,:,:)  ! Hk(nk,norb,norb,nspin,na,na)
                         ! lattice hamiltonian

    double precision, allocatable, public :: &
        dos(:,:,:,:), &  ! dos(nwmodel,norb,nspin,na)
        dosw(:)          ! dos energy grid

    private
contains

    subroutine dmft_lattice_init
        integer :: i,j,ik,nkx,iorb,ispin,ia

        double precision :: dk, kx, ky, ckx, cky

        if (master) then
            write(*,*) "Setting up lattice Hamiltonian..."
        endif

        select case (tbham)
            case (0)
                call read_tb_hamiltonian
            case (1)
                call tbh_model_1
            case (2)
                call tbh_model_2
            case default
                call die("dmft_lattice_init", &
                    "invalid tight-binding hamiltonian type")
        end select
    end subroutine dmft_lattice_init

    subroutine read_tb_hamiltonian
        ! @TODO
        call die("read_tb_hamiltonian", "not implemented yet")

    end subroutine read_tb_hamiltonian

    subroutine tbh_model_1
        ! all degenerate bands on 2D square lattice
        ! hopping = 0.5
        ! bandwidth = 2
        ! half-bandwidth = 1
        integer :: ik,nkx,i,j,iorb,ia,ispin,k
        double precision :: dk, energy, diff, de, kx, ky, &
                            ckx, cky

        allocate(Hk(nk,norb,norb,nspin,na,na))
        hk = cmplx(0.d0,0.d0)

        nkx = int(sqrt(float(nk)))
        dk = 2.0D0*pi/float(nkx-1)

        ik = 0
        do i = 1, nkx
            do j = 1, nkx
                ik = ik + 1
                kx = -pi + dk*float(i-1)
                ky = -pi + dk*float(j-1)
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

        ! density of state
        allocate(dos(nwmodel,norb,nspin,na))
        allocate(dosw(nwmodel))
        de = 3.d0/(nwmodel-1)
        do ia=1,na
            do ispin=1,nspin
                do iorb = 1, norb
                    do i = 1, nwmodel
                        energy = -1.5d0+de*(i-1)
                        dosw(i) = energy
                        dos(i,iorb,ispin,ia) = 0.D0
                        do k=1,nk
                            diff = energy-hk(k,iorb,iorb,ispin,ia,ia)
                            dos(i,iorb,ispin,ia) = dos(i,iorb,ispin,ia) &
                                + exp(-diff*diff/(2.d0*broadening*broadening))&
                                /sqrt(2.d0*pi)/broadening
                        enddo
                        dos(i,iorb,ispin,ia) = dos(i,iorb,ispin,ia)/float(nk)
                    enddo
                enddo
            enddo
        enddo

        call dump_dos
    end subroutine tbh_model_1

    subroutine tbh_model_2
        ! Bethe lattice, infinite coordination
        ! half bandwidth = 1
        integer :: ik,nkx,i,j,iorb,ia,ispin
        double precision :: fac, w, dw, val

        if (na.gt.1) then
            call die("dmft_lattice_init", &
                "Bethe lattice with na>1 is not possible")
        endif

        allocate(dos(nwmodel,norb,nspin,na))
        allocate(dosw(nwmodel))

        fac = 2.d0/pi

        dw = 2.d0/(nwmodel-1)

        do ia=1,na
        do ispin=1,nspin
        do iorb=1,norb
        do i=1,nwmodel
            w = -1.d0 + dw*(i-1)
            dosw(i) = w
            val = 1.d0-w*w
            if (val<0) then
                val = 0.d0
            endif

            dos(i,iorb,ispin,ia) = fac*sqrt(val)
        enddo 
        enddo
        enddo
        enddo

        call dump_dos
    end subroutine tbh_model_2

    subroutine dump_dos
        integer :: ia,ispin,iorb,iw,i
        double precision :: dw,w

        if (master) then
            dw = 2.d0/(nwmodel-1)
            open(unit=IO_LATTICE_DOS, file="dos.dat",status="replace")
            do ia=1,na
            do ispin=1,nspin
            do iorb=1,norb
            do i=1,nwmodel
                write(IO_LATTICE_DOS,"(3I4,2F20.10)") &
                    ia,ispin,iorb,dosw(i),dos(i,iorb,ispin,ia)
            enddo
            enddo
            enddo
            enddo
            close(IO_LATTICE_DOS)
        endif

        call mpi_barrier(comm,mpierr)
    end subroutine dump_dos
end module dmft_lattice
