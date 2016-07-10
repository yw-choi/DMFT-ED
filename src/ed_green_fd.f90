module ed_green_fd
    use dmft_params, only: norb, nspin, na
    use dmft_grid, only: nwloc, omega
    use ed_params, only: nbath, nsite, nsector, &
                         nstep, sectors, nev, eigpair_t
    use ed_basis, only: basis_t, generate_basis
    use lanczos_fd, only: lanczos_iteration_fd
    use ed_operator, only: apply_c
    use ed_hamiltonian, only: generate_hamiltonian
    use numeric_utils, only: mpi_norm

    implicit none

    ! lanczos matrix elements
    double precision, allocatable :: &
        ap(:,:,:,:,:), & ! ap(nstep,nev,norb,nspin,na)
        bp(:,:,:,:,:), & ! bp(nstep,nev,norb,nspin,na)
        an(:,:,:,:,:), & ! an(nstep,nev,norb,nspin,na)
        bn(:,:,:,:,:)    ! bn(nstep,nev,norb,nspin,na)

    double complex, allocatable :: G_cl(:,:,:)

    integer :: ia
contains

    subroutine ed_green_init
        allocate(ap(nstep,nev,norb,nspin,na),bp(nstep,nev,norb,nspin,na))
        allocate(an(nstep,nev,norb,nspin,na),bn(nstep,nev,norb,nspin,na))
        allocate(G_cl(nwloc,norb,nspin))
    end subroutine ed_green_init

    ! calculates the interacting cluster Green's function 
    ! using the Lanczos method.
    subroutine cluster_green_ftn_fd(ia_in,nev_calc,eigpairs)
        integer, intent(in) :: ia_in,nev_calc
        type(eigpair_t), intent(in) :: eigpairs(nev_calc)
        
        ! local variables
        integer :: iev, ispin, iorb, isector, ne_up, ne_down, iw, i
        type(basis_t) :: basis

        ia = ia_in

        G_cl(:,:,:) = cmplx(0.0d0,0.0d0)

        do ispin = 1,nspin
            do iorb = 1,norb
                do iev = 1,nev_calc
                    isector = eigpairs(iev)%sector

                    ne_up = sectors(isector,1)
                    ne_down = sectors(isector,2)

                    call generate_basis(ne_up, ne_down, basis)

                    ! G^+
                    call green_particle(iev,iorb,ispin,basis,eigpairs(iev))

                    ! G^-
                    call green_hole(iev,iorb,ispin,basis,eigpairs(iev))
                
                enddo
            enddo
        enddo

    end subroutine cluster_green_ftn_fd

    subroutine green_particle(iev, iorb, ispin, basis, eigpair)
        type(basis_t), intent(in) :: basis
        type(eigpair_t), intent(in) :: eigpair
        integer, intent(in) :: iev, iorb, ispin

        ! local variables
        double precision, allocatable :: v(:)
        double precision, allocatable :: H(:,:)
        double complex :: gr, z
        type(basis_t) :: basis_out
        integer :: iw

        ! c^+_{iorb,ispin} | eigvec >
        call apply_c( basis, eigpair%vec, 1, iorb, ispin, basis_out, v )

        allocate(H(basis_out%nloc,basis_out%ntot))
        call generate_hamiltonian(basis_out, H)
        call lanczos_iteration_fd( basis_out%nloc, H, v, nstep, &
                                ap(:,iev,iorb,ispin,ia), bp(:,iev,iorb,ispin,ia) )

        do iw = 1,nwloc
            z = cmplx(eigpair%val, omega(iw))
            gr = continued_fraction_p(z, nstep, ap(:,iev,iorb,ispin,ia), &
                                                bp(:,iev,iorb,ispin,ia))
            bp(1,iev,iorb,ispin,ia) = mpi_norm(v,basis_out%nloc)
            gr = gr*bp(1,iev,iorb,ispin,ia)*bp(1,iev,iorb,ispin,ia)*eigpair%prob
            G_cl(iw,iorb,ispin) = G_cl(iw,iorb,ispin) + gr
        enddo

        deallocate(H)
    end subroutine green_particle

    subroutine green_hole(iev, iorb, ispin, basis, eigpair)
        type(basis_t), intent(in) :: basis
        type(eigpair_t), intent(in) :: eigpair
        integer, intent(in) :: iev, iorb, ispin

        ! local variables
        double precision, allocatable :: v(:)
        double precision, allocatable :: H(:,:)
        double precision :: a(nstep), b(nstep)
        double complex :: gr, z
        type(basis_t) :: basis_out
        integer :: iw

        ! c^-_{iorb,ispin} | eigvec >
        call apply_c( basis, eigpair%vec, 2, iorb, ispin, basis_out, v )

        allocate(H(basis_out%nloc,basis_out%ntot))
        call generate_hamiltonian(basis_out, H)
        call lanczos_iteration_fd( basis_out%nloc, H, v, nstep, &
                                an(:,iev,iorb,ispin,ia), bn(:,iev,iorb,ispin,ia) )
        
        do iw = 1,nwloc
            z = cmplx(-eigpair%val, omega(iw))
            gr = continued_fraction_m(z, nstep, an(:,iev,iorb,ispin,ia), &
                                                bn(:,iev,iorb,ispin,ia))
            bn(1,iev,iorb,ispin,ia) = mpi_norm(v,basis_out%nloc)
            gr = gr*bn(1,iev,iorb,ispin,ia)*bn(1,iev,iorb,ispin,ia)*eigpair%prob
            
            g_cl(iw,iorb,ispin) = g_cl(iw,iorb,ispin) + gr
        enddo

        deallocate(H)
    end subroutine green_hole

    ! f =                         1 
    !         ---------------------------------------------
    !                                 b(2)*b(2)  
    !          z - a(1) - ---------------------------------
    !                                          b(3)*b(3)
    !                         z - a(2)  -  ----------------
    !                                       z - a(3) - ... 
    double complex function continued_fraction_p(z,n,a,b) result(f)
        double complex, intent(in) :: z
        integer, intent(in) :: n
        double precision, intent(in) :: a(n), b(n)

        integer :: i 

        f = b(n)*b(n)/(z-a(n))
        do i = n-1, 2, -1
            f = b(i)*b(i)/(z-a(i)-f)
        enddo
        f = 1.d0/(z-a(1)-f)
    end function continued_fraction_p

    ! f =                         1 
    !         ---------------------------------------------
    !                                 b(2)*b(2)  
    !          z + a(1) - ---------------------------------
    !                                          b(3)*b(3)
    !                         z + a(2)  -  ----------------
    !                                       z + a(3) - ... 
    double complex function continued_fraction_m(z,n,a,b) result(f)
        double complex, intent(in) :: z
        integer, intent(in) :: n
        double precision, intent(in) :: a(n), b(n)

        integer :: i 

        f = b(n)*b(n)/(z+a(n))
        do i = n-1, 2, -1
            f = b(i)*b(i)/(z+a(i)-f)
        enddo
        f = 1.d0/(z+a(1)-f)
    end function continued_fraction_m
end module ed_green_fd
