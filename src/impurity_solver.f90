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
    use ed_solver

    implicit none

    character(len=100) :: solver

contains

    subroutine solver_init

        solver = fdf_get("DMFT.Solver", "ED")

        select case(solver)
            case ("ED")
                call ed_init
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select
    end subroutine solver_init

    subroutine solve(iloop,G0,Sigma)
        integer, intent(in) :: iloop
        double complex, intent(in) :: G0(nwloc,norb,nspin)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)

        select case(solver)
            case ("ED")
                call ed_solve(iloop,G0,Sigma)
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select

    end subroutine solve

end module impurity_solver
