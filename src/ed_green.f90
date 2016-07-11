module ed_green
    use dmft_params, only: norb, nspin, na
    use dmft_grid, only: nwloc, omega
    use ed_params, only: nbath, nsite, nsector, &
                         nstep, sectors, nev, eigpair_t
    use ed_basis, only: basis_t, generate_basis
    use lanczos, only: lanczos_iteration
    use ed_operator, only: apply_c
    use numeric_utils, only: mpi_norm, continued_fraction_m, &
                             continued_fraction_p

    implicit none

    type :: g_coeff_t
        integer :: nstep
        integer :: nev

        double precision, allocatable :: ap(:,:)
        double precision, allocatable :: bp(:,:)
        double precision, allocatable :: an(:,:)
        double precision, allocatable :: bn(:,:)

    end type g_coeff_t

    type(g_coeff_t), allocatable :: g_coeffs(:,:,:)
    
    double complex, allocatable :: G_cl(:,:,:)

    integer :: ia
contains

    subroutine ed_green_init
        integer :: ia, ispin, iorb
        allocate(g_coeffs(norb,nspin,na))
        allocate(G_cl(nwloc,norb,nspin))

        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    allocate(g_coeffs(iorb,ispin,ia)%ap(nstep, nev))
                    allocate(g_coeffs(iorb,ispin,ia)%bp(nstep, nev))
                    allocate(g_coeffs(iorb,ispin,ia)%an(nstep, nev))
                    allocate(g_coeffs(iorb,ispin,ia)%bn(nstep, nev))
                enddo
            enddo
        enddo
    end subroutine ed_green_init

    ! calculates the interacting cluster Green's function 
    ! using the Lanczos method.
    subroutine cluster_green_ftn(ia_in,nev_calc,eigpairs)
        integer, intent(in) :: ia_in,nev_calc
        type(eigpair_t), intent(in) :: eigpairs(nev_calc)
        
        ! local variables
        integer :: iev, ispin, iorb, isector, ne_up, ne_down, iw, i
        type(basis_t) :: basis, basis_out
        double precision, allocatable :: v(:)
        double precision :: ap(nstep,nev_calc), bp(nstep,nev_calc), &
                            an(nstep,nev_calc), bn(nstep,nev_calc)
        double complex :: z, gr

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
                    ! 1. c^+|iev>
                    call apply_c(basis, eigpairs(iev)%vec, 1, &
                        iorb, ispin, basis_out, v)

                    ! 2. lanczos coefficients
                    call lanczos_iteration(ia, basis_out, v, nstep, &
                                                ap(:,iev), bp(:,iev))

                    bp(1,iev) = mpi_norm(v, basis_out%nloc)

                    ! 3. green function as a continued fraction
                    do iw=1,nwloc
                        z = cmplx(eigpairs(iev)%val, omega(iw))
                        gr = continued_fraction_p(z, nstep, ap(:,iev), &
                                                            bp(:,iev))
                        gr = gr*eigpairs(iev)%prob
                        G_cl(iw,iorb,ispin) = G_cl(iw,iorb,ispin) + gr
                    enddo

                    ! G^-
                    ! 1. c^-_{iorb,ispin} | eigvec >
                    call apply_c( basis, eigpairs(iev)%vec, 2, &
                        iorb, ispin, basis_out, v )

                    ! 2. lanczos coefficients
                    call lanczos_iteration(ia, basis_out, v, nstep, &
                                               an(:,iev), bn(:,iev))
                    bn(1,iev) = mpi_norm(v, basis_out%nloc)

                    ! 3. green function as a continued fraction
                    do iw = 1,nwloc
                        z = cmplx(-eigpairs(iev)%val, omega(iw))
                        gr = continued_fraction_m(z, nstep, an(:,iev), &
                                                            bn(:,iev))
                        gr = gr*eigpairs(iev)%prob
                        
                        g_cl(iw,iorb,ispin) = g_cl(iw,iorb,ispin) + gr
                    enddo
                enddo

                g_coeffs(iorb,ispin,ia)%nstep = nstep
                g_coeffs(iorb,ispin,ia)%nev   = nev_calc
                g_coeffs(iorb,ispin,ia)%ap    = ap
                g_coeffs(iorb,ispin,ia)%bp    = bp
                g_coeffs(iorb,ispin,ia)%an    = an
                g_coeffs(iorb,ispin,ia)%bn    = bn

            enddo
        enddo

    end subroutine cluster_green_ftn

end module ed_green
