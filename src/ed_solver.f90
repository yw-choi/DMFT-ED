module ed_solver
    use mpi
    use dmft_params, only: norb
    use matsubara_grid, only: nwloc
    use ed_projection, only: project_to_impurity_model
    use ed_params, only: ek,vk

    implicit none


contains

    subroutine ed_solve(G0,Sigma)
        double complex, intent(in) :: G0(norb,nwloc)
        double complex, intent(out) :: Sigma(norb,nwloc)
        double complex :: G0_cl(norb,nwloc) ! cluster-projected Weiss field

        integer :: i,j

        call project_to_impurity_model(G0,ek,vk)



    end subroutine ed_solve

end module ed_solver
