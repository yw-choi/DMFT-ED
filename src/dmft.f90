module dmft
!==============================================================================
! controller module for dmft calculation.
! this module takes care of the flow of the program,
! and the actual works are done by other modules called by this module.
!==============================================================================
    use dmft_grid
    use dmft_params
    use dmft_lattice
    use dmft_green, only: G, G0, Sigma, G_prev, &
                          dmft_green_init, &
                          update_local_green_ftn, &
                          update_weiss_ftn, &
                          dump_green

    use impurity_solver, only: impurity_solver_init, solve, solver_post_processing

    use utils
    use mpi
    use timer, only: t1_loop, t2_loop, print_elapsed_time, t1_run
    implicit none

    integer ::  &
        iloop     ! DMFT loop count

    logical :: &
        converged  ! is DMFT loop converged?


contains

    subroutine dmft_init
        ! Read general DMFT parameters independent of solver.
        call dmft_params_read
        call dmft_grid_init
        call dmft_lattice_init
        call dmft_green_init

        call impurity_solver_init
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

            call update_local_green_ftn
            call update_weiss_ftn
            call dump_green 

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

    subroutine loop_end
        integer :: iw, ia, iorb
        double precision :: diff, diffsum

        ! Test DMFT convergence.
        diff = sum(abs(G_prev(:,:,:,:)-G(:,:,:,:)))
        diff = diff / (na*norb*nspin*nwloc)

        call mpi_allreduce(diff,                 &
                           diffsum,              &
                           1,                    &
                           mpi_double_precision, &
                           mpi_sum,              &
                           comm,                 &
                           mpierr)

        if (diffsum < scf_tol) then
            converged = .true.
        endif

        if (master) then
            write(*,*)
            write(*,"(1x,a,I4,a,E12.5)") "loop ",iloop, &
                                         " done. scf diff = ",diffsum
            call print_elapsed_time(" Elapsed time in a loop     ",&
                                    t1_loop,t2_loop)
            call print_elapsed_time(" Elapsed time from the start",&
                                    t1_run,t2_loop)
            if (converged) then
                write(*,*) "DMFT loop has converged within ",iloop," iterations."
            endif
            write(*,*)
        endif

    end subroutine loop_end
end module dmft
