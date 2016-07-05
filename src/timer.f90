module timer

    implicit none

    double precision :: &
        t1_run, t2_run,                 &  ! Total Running Time
        t1_solver, t2_solver,           &  ! solver
        t1_update_loop, t2_update_loop, &  ! updating lattice green/weiss field
        t1_diag_loop, t2_diag_loop,     &  ! Diagonalization time in a loop
        t1_green_loop, t2_green_loop       ! Green's ftn calculation 

contains

    double precision function minutes(t1,t2)
        double precision :: t1,t2
        minutes = (t2-t1)/60.d0
    end function minutes

    subroutine print_elapsed_time(msg, t1, t2)
        character(len=*) :: msg
        double precision :: t1, t2

        write(*,"(2a,3x,f10.5,3x,a)") msg, " : " , minutes(t1,t2), " min."

    end subroutine print_elapsed_time

end module timer
