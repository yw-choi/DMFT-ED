module dmft
    use matsubara_grid
    implicit none

    integer ::  &
        norb,   & ! number of impurity orbitals
        nloop,  & ! maximum number of DMFT loop
        nk        ! number of kpoints in band structure (@TODO to be refactored out)

    double precision :: &
        U,            & ! U
        Jex,          & ! J
        Mu,           & ! chemical potential
        beta,         & ! inverse temperature
        broadening,   & ! inverse temperature
        scf_tol         ! DMFT SCF tolerance

    double complex, allocatable :: &
        G_prev(:,:),  & ! G_prev(norb,nwloc) local Green's function of the previous step
        G(:,:),       & ! G(norb,nwloc)      local Green's function of the current step
        G0(:,:),      & ! G0(norb,nwloc)     Weiss field of the current step
        Sigma(:,:)      ! Sigma(norb,nwloc)  Self-energy of the current step

contains

    ! Read general DMFT parameters independent of solver.
    subroutine read_dmft_options
        use fdf
        use mpi
        use utils
        use matsubara_grid, only: nw
        
        implicit none

        character(len=200) :: msg

        if (master) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,'(4x,a)') "DMFT Parameters"
            write(6,'(a)') repeat("=",80)
        endif

        norb = fdf_get("DMFT.Norb",0)
        if (norb.eq.0.and.master) then
            call die("read_dmft_options", "Number of orbitals must be larger than 0.")
        endif

        U = fdf_get("DMFT.U", 1.0D0)
        Jex = fdf_get("DMFT.J",0.3D0) 
        Mu = fdf_get("DMFT.Mu",1.0D0)
        beta = fdf_get("DMFT.beta", 100.0D0)

        Nloop = fdf_get("DMFT.MaxIterations", 100)
        scf_tol = fdf_get("DMFT.ScfTolerance", 1d-4)

        Nw = fdf_get("DMFT.Matsubara.Nw",2000)

        nk = fdf_get("DMFT.Nk",10000)
        broadening = fdf_get("DMFT.SpectralBroadening", 0.02D0)

        if (master) then
            write(msg,*) "Number of impurity orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', Norb
            write(msg,*) "U"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', U
            write(msg,*) "J"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', Jex
            write(msg,*) "chemical potential"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', Mu
            write(msg,*) "beta (inverse temperature)"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', beta
            write(msg,*) "Maximum DMFT iterations"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', Nloop
            write(msg,*) "DMFT SCF tolerance"
            write(6,'(3x,a40,2x,a,2x,ES8.1)') msg, '=', scf_tol
            write(msg,*) "Number of Matsubara frequencies"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nw
            write(msg,*) "Number of k-points in band structure"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nk
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif
    end subroutine read_dmft_options
end module dmft
