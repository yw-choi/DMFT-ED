module dmft_params
! =============================================================================
! DMFT parameters that do not change during the run.
! =============================================================================

    integer ::      &
        nspin,      & ! number of spin components 
        na,         & ! number of atoms in unit cell
        norb,       & ! number of impurity orbitals
        nloop,      & ! maximum number of DMFT loop
        nw,         & ! Total number of matsubara frequencies
        nwreal,     & ! Total number of real frequencies
        nwmodel,    & ! # of frequency point used in model dos calculation
        nk,         & ! number of k-points in lattice hamiltonian
        tbham         ! type of tight-binding hamiltonian
                      ! 0 : read from file
                      ! 1 : 2d square lattice
                      ! 2 : bethe lattice, infinite dimension (circular dos)

    double precision :: &
        U,            & ! U
        Up,           & ! U'
        Jex,          & ! J
        Jp,           & ! J'
        Mu,           & ! chemical potential
        beta,         & ! inverse temperature
        broadening,   & ! inverse temperature
        scf_tol,      & ! DMFT SCF tolerance
        wrmin,        & ! real frequency min
        wrmax           ! real frequency max

    logical :: &
        read_gf_from_file

contains

    ! Read general DMFT parameters independent of solver.
    subroutine dmft_params_read
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
            call die("dmft_params_read", "Number of atoms must be larger than 0.")
        endif

        if (nspin.gt.2.and.master) then
            call die("dmft_params_read", "currently, nspin=1,2 is implemented.")
        endif

        norb = fdf_get("DMFT.Norb",0)
        if (norb.eq.0.and.master) then
            call die("dmft_params_read", "Number of orbitals must be larger than 0.")
        endif

        U = fdf_get("DMFT.U", 1.0D0)
        Jex = fdf_get("DMFT.J",0.3D0) 
        Up = fdf_get("DMFT.Up",-1.0d0)
        Jp = fdf_get("DMFT.Jp",-1.0d0) 

        if (Up<0) then
            Up = U - 2*Jex
        endif

        if (Jp<0) then
            Jp = Jex
        endif

        Mu = fdf_get("DMFT.Mu",1.0D0)
        beta = fdf_get("DMFT.beta", 100.0D0)

        Nloop = fdf_get("DMFT.MaxIterations", 100)
        scf_tol = fdf_get("DMFT.ScfTolerance", 1d-4)

        nw = fdf_get("DMFT.Matsubara.Nw",1000)
        nwreal = fdf_get("DMFT.Real.Nw",4000)
        wrmin = fdf_get("DMFT.Real.wmin",-3.d0)
        wrmax = fdf_get("DMFT.Real.wmax", 3.d0)
        nwmodel = fdf_get("DMFT.ModelDOSNw",2000)
        nk = fdf_get("DMFT.Nk",40000)
        broadening = fdf_get("DMFT.SpectralBroadening", 0.02D0)

        tbham = fdf_get("DMFT.TBHamiltonian", 1)

        read_gf_from_file = fdf_get("DMFT.GreenFtnFromFile", .false.)

        if (master) then
            write(msg,*) "Number of atoms in a unit cell"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', na
            write(msg,*) "Number of spin components"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', nspin
            write(msg,*) "Number of impurity orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', Norb
            write(msg,*) "U"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', U
            write(msg,*) "Up"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', Up
            write(msg,*) "J"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', Jex
            write(msg,*) "Jp"
            write(6,'(3x,a40,2x,a,2x,F8.3)') msg, '=', Jp
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
            write(msg,*) "Tight-binding Hamiltonian"
            write(6,'(3x,a40,2x,a,2x,I8)') msg, '=', tbham
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif

    end subroutine dmft_params_read

end module dmft_params
