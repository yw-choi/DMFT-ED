module impurity_solver
!==============================================================================
! Wrapper module for impurity solvers. 
! This module separates the specific implementation of the impurity solver with
! DMFT module. 
!==============================================================================
    use ed_solver
    use fdf
    use utils
    implicit none

    character(len=100) :: solver

contains

    subroutine solver_init

        solver = fdf_get("DMFT.Solver", "ED")

        select case(solver)
            case ("ED")
                call ed_read_options
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select
    end subroutine solver_init

    subroutine solve

        select case(solver)
            case ("ED")
                call ed_solve
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select

    end subroutine solve

end module impurity_solver
