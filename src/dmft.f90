module dmft
    use matsubara_grid
    use impurity_solver, only: solver_init, solve, solver_post_processing
    use dmft_params
    use dmft_lattice
    use utils
    use mpi
    use timer, only: t1_loop, t2_loop, print_elapsed_time, t1_run
    implicit none

    integer ::  &
        iloop     ! DMFT loop count

    logical :: &
        converged  ! is DMFT loop converged?

    double complex, allocatable :: &
        G_prev(:,:,:,:),   & ! G_prev(nwloc,norb,nspin,na) local Green's function of the previous step
        G(:,:,:,:),        & ! G(nwloc,norb,nspin,na)      local Green's function of the current step
        G0(:,:,:,:),       & ! G0(nwloc,norb,nspin,na)     Weiss field of the current step
        Sigma(:,:,:,:)       ! Sigma(nwloc,norb,nspin,na)  Self-energy of the current step

contains

    subroutine dmft_init
        ! Read general DMFT parameters independent of solver.
        call read_dmft_params

        call setup_matsubara_grid

        ! impurity-solver-specific initialization
        call solver_init

        allocate(G_prev(nwloc,norb,nspin,na)) 
        allocate(G(nwloc,norb,nspin,na))
        allocate(G0(nwloc,norb,nspin,na))
        allocate(Sigma(nwloc,norb,nspin,na))

        G_prev = cmplx(0.0d0,0.0d0)
        G      = cmplx(0.0d0,0.0d0)
        G0     = cmplx(0.0d0,0.0d0)
        Sigma  = cmplx(0.0d0,0.0d0)

        call setup_lattice_hamiltonian

        ! Setting up the initial Weiss field
        call update_local_green_ftn
        call update_weiss_ftn

    end subroutine dmft_init

    subroutine dmft_loop
        integer :: ispin, ia

        if (master) then
            write(6,"(a)") repeat("=",80)
            write(*,*)
            write(*,*) "Start of DMFT SCF loop"
            write(*,*)
            write(6,"(a)") repeat("=",80)
        endif

        converged = .false.
        iloop = 0
        do while(.not.converged.and.iloop.le.nloop)
            iloop = iloop + 1
            if (master) then
                write(6,"(a)") repeat("=",80)
                write(*,*)
                write(*,*) "DMFT Loop ", iloop
                write(*,*)
                write(6,"(a)") repeat("=",80)
            endif

            t1_loop = mpi_wtime(mpierr)

            do ia=1,na
                if (master) then
                    write(*,"(1x,a,I1,a,I2,a)") "impurity problem for ia = ",ia
                endif

                call solve(iloop,ia,G0(:,:,:,ia), Sigma(:,:,:,ia))
            enddo

            if (dump_loop_data) call dump_data

            call update_local_green_ftn
            call update_weiss_ftn

            t2_loop = mpi_wtime(mpierr)
            call loop_end
        enddo

        if (master) then
            write(*,"(a)") repeat("=",80)
            write(*,*)
            write(*,*) "End of DMFT SCF loop"
            write(*,*)
            write(*,"(a)") repeat("=",80)
        endif
    end subroutine dmft_loop

    subroutine dmft_post_processing
        integer :: iw, iorb, ispin, ia
        character(len=100) fn
        double precision :: zq

        if (master) then
            ! dump final G0 G Sigma (only low frequency data)
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        write(fn,"(a,I1.1,a1,I1.1,a1,I1.1,a)") &
                            "green_ftn.",ia,".",ispin,".",iorb,".dat"
                        open(unit=11,file=fn,form="formatted",status="replace")
                        do iw=1,nwloc
                            write(11,"(3F15.6)") omega(iw),real(G(iw,iorb,ispin,ia)),aimag(G(iw,iorb,ispin,ia))
                        enddo
                        close(11)
                        write(fn,"(a,I1.1,a1,I1.1,a1,I1.1,a)") &
                            "self_energy.",ia,".",ispin,".",iorb,".dat"
                        open(unit=11,file=fn,form="formatted",status="replace")
                        do iw=1,nwloc
                            write(11,"(3F15.6)") omega(iw),real(Sigma(iw,iorb,ispin,ia)),aimag(Sigma(iw,iorb,ispin,ia))
                        enddo
                        close(11)
                        write(fn,"(a,I1.1,a1,I1.1,a1,I1.1,a)") &
                            "weiss_field.",ia,".",ispin,".",iorb,".dat"
                        open(unit=11,file=fn,form="formatted",status="replace")
                        do iw=1,nwloc
                            write(11,"(3F15.6)") omega(iw),real(G0(iw,iorb,ispin,ia)),aimag(G0(iw,iorb,ispin,ia))
                        enddo
                        close(11)
                    enddo
                enddo
            enddo

            write(*,*) 
            write(*,*) "Quasiparticle weight"
            write(*,*) 
            write(*,*) "    Z(ia,ispin,iorb)"
            write(*,*) "    ---------------------------"
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        zq = aimag(sigma(1,iorb,ispin,ia))&
                            /omega(1)                        
                        zq = 1.d0/(1.d0-zq)
                        write(*,"(a,I2,a,I4,a,I3,a,F12.6)") &
                            "     Z(",ia,",",ispin,",",iorb,",) = ",zq
                    enddo
                enddo
            enddo
            write(*,*) 
        endif

        call mpi_barrier(comm,mpierr)

        call solver_post_processing

    end subroutine dmft_post_processing

    subroutine dmft_finalize
        deallocate(G_prev,G,G0,Sigma)
    end subroutine dmft_finalize

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

    ! Test DMFT convergence.
    subroutine loop_end
        integer :: iw, ia, iorb
        double precision :: diff, diffsum

        diff = sum(abs(G_prev(:,:,:,:)-G(:,:,:,:)))
        diff = diff / (na*norb*nspin*nwloc)

        call mpi_allreduce(diff,diffsum,1,mpi_double_precision,mpi_sum,comm,mpierr)

        if (master) then
            write(*,*)
            write(*,"(1x,a,I4,a,E12.5)") "loop ",iloop," done. scf diff = ",diffsum
            call print_elapsed_time(" Elapsed time in a loop     ",&
                                    t1_loop,t2_loop)
            call print_elapsed_time(" Elapsed time from the start",&
                                    t1_run,t2_loop)
            write(*,*)
        endif

        if (diffsum < scf_tol) then
            converged = .true.
        endif
    end subroutine loop_end

    subroutine dump_data
    ! dump current green ftn, weiss ftn, self energy
        use utils
        use io_units
        integer :: iw, iorb, ia, ispin
        character(len=100) :: fn
        double precision :: zq

        if (master) then
            write(fn,"(a5,I3.3,a4)") "dmft.",iloop,".dat"
            open(unit=IO_DEBUG_DUMP,file=fn,status="replace")
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        do iw=1,nwloc
                            write(IO_DEBUG_DUMP,"(7F16.8)") omega(iw), &
                                real(G(iw,iorb,ispin,ia)), aimag(G(iw,iorb,ispin,ia)), &
                                real(G0(iw,iorb,ispin,ia)), aimag(G0(iw,iorb,ispin,ia)), &
                                real(Sigma(iw,iorb,ispin,ia)), aimag(Sigma(iw,iorb,ispin,ia))
                        enddo
                        write(IO_DEBUG_DUMP,*)
                    enddo
                enddo
            enddo
            close(IO_DEBUG_DUMP)
        endif
        
        call mpi_barrier(comm,mpierr)
    end subroutine dump_data
end module dmft
