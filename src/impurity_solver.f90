module impurity_solver
    use ed_solver
    use fdf
    use utils
    implicit none

    character(len=100) :: solver

contains

    subroutine solve

        select case(solver)
            case ("ED")
                call ed_solve
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select

    end subroutine solve

    subroutine read_solver_options

        solver = fdf_get("DMFT.Solver", "ED")

        select case(solver)
            case ("ED")
                call ed_read_options
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select
    end subroutine read_solver_options

end module impurity_solver
