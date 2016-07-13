module ed_green
    use dmft_params, only: norb, nspin, na
    use dmft_grid, only: nwloc, omega
    use ed_params, only: nbath, nsite, nsector, &
                         nstep, sectors, nev
    use ed_basis, only: basis_t, generate_basis, dealloc_basis
    use lanczos, only: lanczos_iteration
    use ed_operator, only: apply_c
    use numeric_utils, only: mpi_norm, continued_fraction_m, &
                             continued_fraction_p
    use ed_eigpair, only: eigpair_t, nev_calc, eigpairs
    use alloc, only: re_alloc, de_alloc

    use mpi

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

    double precision, pointer :: vec_out(:), vec_all(:)
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
    subroutine cluster_green_ftn(ia)
        integer, intent(in) :: ia
        
        ! local variables
        integer :: iev, ispin, iorb, isector, ne_up, ne_down, iw, i
        type(basis_t) :: basis
        double complex, allocatable :: gdiag_up(:), gdiag_down(:)

        allocate(gdiag_up(nwloc),gdiag_down(nwloc))
        nullify(vec_all, vec_out)

        G_cl(:,:,:) = cmplx(0.0d0,0.0d0)

        do iorb = 1,norb
            do iev = 1,nev_calc
                ispin = 1
                if (master) then
                    write(*,"(1x,A,A,I2,A,I2,A,I2,A,I2,A)") &
                        "Cluster Green's function for (ia,ispin,iorb,iev)=",&
                        "(",ia,",",ispin,",",iorb,",",iev,")"

                endif
                isector = eigpairs(iev)%sector

                ne_up   = sectors(isector,1)
                ne_down = sectors(isector,2)

                call generate_basis(ne_up, ne_down, basis)

                call re_alloc(vec_all, 1, basis%ntot, 'ed_green', 'vec_all')

                call green_diag(ia,iorb,iev,ispin,basis,gdiag_up)

                ! in case of paramagnetic case with same up/down spin
                ! there is no need to calculate G_down
                if (nspin==2.or.ne_up/=ne_down) then
                    ispin = 2
                    if (master) then
                        write(*,"(1x,A,A,I2,A,I2,A,I2,A,I2,A)") &
                            "Cluster Green's function for (ia,ispin,iorb,iev)=",&
                            "(",ia,",",ispin,",",iorb,",",iev,")"

                    endif
                    call green_diag(ia,iorb,iev,ispin,basis,gdiag_down)
                endif

                if (nspin==1.and.ne_up==ne_down) then
                    G_cl(:,iorb,1) = G_cl(:,iorb,1) + gdiag_up(:)
                else if (nspin==1.and.ne_up/=ne_down) then
                    G_cl(:,iorb,1) = G_cl(:,iorb,1) + &
                        0.5d0*(gdiag_up(:)+gdiag_down(:))
                else
                    G_cl(:,iorb,1) = G_cl(:,iorb,1) + gdiag_up(:)
                    G_cl(:,iorb,2) = G_cl(:,iorb,2) + gdiag_down(:)
                endif
            enddo
        enddo

        call dealloc_basis(basis)
        call de_alloc(vec_all, 'ed_green', 'vec_all')
        call de_alloc(vec_out, 'ed_green', 'vec_out')

    end subroutine cluster_green_ftn

    subroutine green_diag(ia,iorb,iev,ispin,basis,gdiag)
        integer, intent(in) :: ia,iorb,iev,ispin
        type(basis_t), intent(in) :: basis
        double complex, intent(out) :: gdiag(nwloc)

        integer :: iw
        type(basis_t) :: basis_out
        double precision, allocatable :: ap(:,:), bp(:,:), &
                                         an(:,:), bn(:,:)
        double complex :: z, gr

        allocate(ap(nstep,nev_calc),bp(nstep,nev_calc))
        allocate(an(nstep,nev_calc),bn(nstep,nev_calc))

        gdiag = cmplx(0.d0,0.d0)

        ! G^+
        ! 1. c^+|iev>
        if (ispin == 1) then
            call generate_basis( basis%ne_up+1, basis%ne_down, &
                                 basis_out)
        else
            call generate_basis( basis%ne_up, basis%ne_down+1, &
                                 basis_out)
        endif

        call re_alloc(vec_out, 1, basis_out%nloc, 'ed_green', 'vec_out')

        call apply_c(basis, eigpairs(iev)%vec, vec_all, 1, &
            iorb, ispin, basis_out, vec_out)

        ! 2. lanczos coefficients
        call lanczos_iteration(ia, basis_out, vec_out, nstep, &
                                    ap(:,iev), bp(:,iev))

        bp(1,iev) = mpi_norm(vec_out, basis_out%nloc)

        ! 3. green function as a continued fraction
        do iw=1,nwloc
            z = cmplx(eigpairs(iev)%val, omega(iw))
            gr = continued_fraction_p(z, nstep, ap(:,iev), &
                                                bp(:,iev))
            gr = gr*eigpairs(iev)%prob
            gdiag(iw) = gdiag(iw)+gr
        enddo

        ! G^-
        ! 1. c^-_{iorb,ispin} | eigvec >
        if (ispin == 1) then
            call generate_basis( basis%ne_up-1, basis%ne_down, &
                                 basis_out)
        else
            call generate_basis( basis%ne_up, basis%ne_down-1, &
                                 basis_out)
        endif

        call re_alloc(vec_out, 1, basis_out%nloc, 'ed_green', 'vec_out')

        call apply_c( basis, eigpairs(iev)%vec, vec_all, 2, &
            iorb, ispin, basis_out, vec_out )

        ! 2. lanczos coefficients
        call lanczos_iteration(ia, basis_out, vec_out, nstep, &
                                   an(:,iev), bn(:,iev))
        bn(1,iev) = mpi_norm(vec_out, basis_out%nloc)

        ! 3. green function as a continued fraction
        do iw = 1,nwloc
            z = cmplx(-eigpairs(iev)%val, omega(iw))
            gr = continued_fraction_m(z, nstep, an(:,iev), &
                                                bn(:,iev))
            gr = gr*eigpairs(iev)%prob
            
            gdiag(iw) = gdiag(iw)+gr
        enddo

        g_coeffs(iorb,ispin,ia)%nstep = nstep
        g_coeffs(iorb,ispin,ia)%nev   = nev_calc

        g_coeffs(iorb,ispin,ia)%ap(:,1:nev_calc) = ap(:,1:nev_calc)
        g_coeffs(iorb,ispin,ia)%bp(:,1:nev_calc) = bp(:,1:nev_calc)

        g_coeffs(iorb,ispin,ia)%an(:,1:nev_calc) = an(:,1:nev_calc)
        g_coeffs(iorb,ispin,ia)%bn(:,1:nev_calc) = bn(:,1:nev_calc)

        deallocate(ap,bp,an,bn)
        call dealloc_basis(basis_out)
    end subroutine green_diag
end module ed_green
