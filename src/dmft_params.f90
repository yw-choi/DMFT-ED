module dmft_params
! =============================================================================
! DMFT parameters that do not change during the run.
! =============================================================================

    integer ::  &
        nspin,  & ! number of spin components 
        na,    & ! number of atoms in unit cell
        norb,   & ! number of impurity orbitals
        nloop,  & ! maximum number of DMFT loop
        nw,     & ! Total number of matsubara frequencies
        nwreal, & ! Total number of real frequencies
        nk

    double precision :: &
        U,            & ! U
        Jex,          & ! J
        Mu,           & ! chemical potential
        beta,         & ! inverse temperature
        broadening,   & ! inverse temperature
        scf_tol         ! DMFT SCF tolerance

contains

    ! Read general DMFT parameters independent of solver.
    subroutine read_dmft_params
        use fdf
        use mpi
        use utils

        character(len=200) :: msg

        if (master) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,'(4x,a)') "DMFT Parameters"
            write(6,'(a)') repeat("=",80)
        endif

        na  = fdf_get("DMFT.Na",0)
        nspin = fdf_get("DMFT.nspin",1)
        if (na.eq.0.and.master) then
            call die("read_dmft_params", "Number of atoms must be larger than 0.")
        endif

        if (nspin.gt.2.and.master) then
            call die("read_dmft_params", "currently, nspin=1,2 is implemented.")
        endif

        norb = fdf_get("DMFT.Norb",0)
        if (norb.eq.0.and.master) then
            call die("read_dmft_params", "Number of orbitals must be larger than 0.")
        endif

        U = fdf_get("DMFT.U", 1.0D0)
        Jex = fdf_get("DMFT.J",0.3D0) 
        Mu = fdf_get("DMFT.Mu",1.0D0)
        beta = fdf_get("DMFT.beta", 100.0D0)

        Nloop = fdf_get("DMFT.MaxIterations", 100)
        scf_tol = fdf_get("DMFT.ScfTolerance", 1d-4)

        nw = fdf_get("DMFT.Matsubara.Nw",2000)
        nwreal = fdf_get("DMFT.Real.Nw",4000)
        nk = fdf_get("DMFT.Nk",10000)
        broadening = fdf_get("DMFT.SpectralBroadening", 0.02D0)

        if (master) then
            write(msg,*) "Number of atoms in unit cell"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', na
            write(msg,*) "Number of spin components"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nspin
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
            write(msg,*) "Number of k-points in band structure"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nk
            write(msg,*) "Number of Matsubara Frequencies"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nw
            write(msg,*) "Number of Real Frequencies"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nwreal
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif

    end subroutine read_dmft_params

end module dmft_params
