module impurity_solver
!==============================================================================
! Wrapper module for impurity solvers. 
! This module separates the specific implementation of the impurity solver with
! DMFT module. 
!==============================================================================
    use fdf
    use utils
    use dmft_params
    use matsubara_grid
    implicit none

    character(len=100) :: solver

contains

    subroutine solver_init
        use ed_params

        solver = fdf_get("DMFT.Solver", "ED")

        select case(solver)
            case ("ED")
                call ed_read_options
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select
    end subroutine solver_init

    subroutine solve(G0,Sigma)
        use ed_solver
        double complex :: G0(norb,nwloc), Sigma(norb,nwloc)

        select case(solver)
            case ("ED")
                call ed_solve(G0,Sigma)
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select

    end subroutine solve

end module impurity_solver
