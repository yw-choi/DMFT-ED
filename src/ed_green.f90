module ed_green
    use dmft_params, only: norb, nspin
    use matsubara_grid, only: nwloc, omega
    use ed_params, only: nbath, nsite, nsector, &
                         nstep, sectors, eigpair_t, nev
    use ed_basis, only: basis_t, generate_basis
    use lanczos, only: lanczos_iteration
    use ed_operator, only: apply_c
    use ed_hamiltonian, only: generate_hamiltonian
    use numeric_utils, only: mpi_norm

    implicit none

    ! lanczos matrix elements
    double precision, allocatable :: &
        ap(:,:,:,:), & ! ap(nstep,nev,norb,nspin)
        bp(:,:,:,:), & ! bp(nstep,nev,norb,nspin)
        an(:,:,:,:), & ! an(nstep,nev,norb,nspin)
        bn(:,:,:,:)    ! bn(nstep,nev,norb,nspin)

    double complex, allocatable :: G_cl(:,:,:)
    double precision, allocatable :: H(:,:)
contains

    subroutine ed_green_init
        allocate(ap(nstep,nev,norb,nspin),bp(nstep,nev,norb,nspin))
        allocate(an(nstep,nev,norb,nspin),bn(nstep,nev,norb,nspin))
        allocate(G_cl(nwloc,norb,nspin))
    end subroutine ed_green_init

    ! calculates the interacting cluster Green's function 
    ! using the Lanczos method.
    subroutine cluster_green_ftn(nev_calc,eigpairs)
        integer, intent(in) :: nev_calc
        type(eigpair_t), intent(in) :: eigpairs(nev_calc)
        
        ! local variables
        integer :: iev, ispin, iorb, isector, ne_up, ne_down, iw, i
        type(basis_t) :: basis

        G_cl(:,:,:) = cmplx(0.0d0,0.0d0)

        do ispin = 1,nspin
            do iorb = 1,norb
                ! G(iorb,ispin) = 
                !          sum_{iev} G^+(iev,iorb,ispin,iev) + G^-(iev,iorb,ispin)
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

        ! open(unit=91,file="gcl.dat",form="formatted")
        ! do iw = 1,nwloc
        !     write(91,"(3F20.8)") omega(iw), real(g_cl(iw,1,1)), aimag(g_cl(iw,1,1))
        ! enddo
        ! close(91)
        ! open(unit=91,file="apbp.dat",form="formatted")
        ! do iev=1,nev_calc
        !     write(91,*) "ap       bp"
        !     do i=1,nstep
        !         write(91,"(2F12.6)") ap(i,iev,1,1),bp(i,iev,1,1)
        !     enddo
        !     write(91,*)
        !     write(91,*) "an       bn"
        !     do i=1,nstep
        !         write(91,"(2F12.6)") an(i,iev,1,1),bn(i,iev,1,1)
        !     enddo
        !     write(91,*)

        ! enddo
        ! close(91)
        ! stop
    end subroutine cluster_green_ftn

    subroutine multiply_h(n,x,y)
        ! @TODO should be parallelized
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: y(n)

        integer :: i,j

        do i=1,n
            y(i) = 0.0d0
            do j=1,n
                ! note that h(j,i) = h(i,j)
                y(i) = y(i) + h(j,i)*x(j) 
            enddo
        enddo
    end subroutine multiply_h

    subroutine green_particle(iev, iorb, ispin, basis, eigpair)
        type(basis_t), intent(in) :: basis
        type(eigpair_t), intent(in) :: eigpair
        integer, intent(in) :: iev, iorb, ispin

        ! local variables
        double precision, allocatable :: v(:)
        double complex :: gr, z
        type(basis_t) :: basis_out
        integer :: iw

        ! c^+_{iorb,ispin} | eigvec >
        call apply_c( basis, eigpair%vec, 1, iorb, ispin, basis_out, v )

        ! @TODO should be parallelized
        allocate(H(basis_out%nloc,basis_out%ntot))
        call generate_hamiltonian(basis_out, H)
        call lanczos_iteration( multiply_h, basis_out%nloc, v, nstep, &
                                ap(:,iev,iorb,ispin), bp(:,iev,iorb,ispin) )

        do iw = 1,nwloc
            z = cmplx(eigpair%val, omega(iw))
            gr = continued_fraction_p(z, nstep, ap(:,iev,iorb,ispin), &
                                                bp(:,iev,iorb,ispin))
            bp(1,iev,iorb,ispin) = mpi_norm(v,basis_out%nloc)
            gr = gr*bp(1,iev,iorb,ispin)*bp(1,iev,iorb,ispin)*eigpair%prob
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
        double precision :: a(nstep), b(nstep)
        double complex :: gr, z
        type(basis_t) :: basis_out
        integer :: iw

        ! c^-_{iorb,ispin} | eigvec >
        call apply_c( basis, eigpair%vec, 2, iorb, ispin, basis_out, v )

        ! @TODO should be parallelized
        allocate(H(basis_out%nloc,basis_out%ntot))
        call generate_hamiltonian(basis_out, H)
        call lanczos_iteration( multiply_h, basis_out%nloc, v, nstep, &
                                an(:,iev,iorb,ispin), bn(:,iev,iorb,ispin) )
        
        do iw = 1,nwloc
            z = cmplx(-eigpair%val, omega(iw))
            gr = continued_fraction_m(z, nstep, an(:,iev,iorb,ispin), &
                                                bn(:,iev,iorb,ispin))
            bn(1,iev,iorb,ispin) = mpi_norm(v,basis_out%nloc)
            gr = gr*bn(1,iev,iorb,ispin)*bn(1,iev,iorb,ispin)*eigpair%prob
            
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
end module ed_green
