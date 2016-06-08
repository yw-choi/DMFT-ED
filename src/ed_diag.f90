module ed_diag

    use mpi
    use ed_params, only: nev, nsector

    implicit none

    type eigval_t
        double precision :: val   ! eigenvalue 
        double precision :: prob  ! exp(-beta*val)/Z
        integer  :: sector        ! sector index
        integer  :: level         ! level index within the sector
    end type eigval_t

    type(eigval_t), allocatable :: eigval(:)
    integer :: nev_calc

    private
contains

    subroutine diagonalize
        integer :: isector, i, ne_up, ne_down
        double precision, allocatable :: eigvec(:,:)

        double precision :: eigval_sec(nev), eigval_all(nev*nsector)
        double precision :: pev(nev*nsector), Z
        integer :: ind(nsector*nev)
        type(basis_t) :: basis
        character(len=100) :: timerstr

        do isector=1,nsector
            ne_up = nup(isector)
            ne_down = nelec(isector) - nup(isector)

            if (node.eq.0) then
                write(*,"(x,a,I3)") "ed_solve: sector ",isector
            endif

            basis = generate_basis( ne_up, ne_down)

            allocate(eigvec(basis%nloc,nev))

            call timestamp2("diag start")
            write(timerstr,"(a6,I1)") "sector",isector
            call timer(timerstr,1)
            call diag(basis,eigval_sec,eigvec)
            call timer(timerstr,2)
            call timestamp2("diag end  ")

            do i=1,nev
                eigval_all((isector-1)*nev+i) = eigval_sec(i)
                call export_eigvec(isector,i,node,basis%nloc,nev,eigvec(:,i))
            enddo

            deallocate(eigvec)
        enddo

        call sort(nev*nsector,eigval_all,ind)
        call boltzmann_factor(eigval_all,nev*nsector,beta,pev)

        nev_calc = 0

        do i = 1, nev
            if(pev(ind(i)).gt.PROB_THRESHOLD) nev_calc = nev_calc + 1
        enddo

        if (allocated(eigval)) deallocate(eigval)

        allocate(eigval(nev_calc))
        Z = 0.0D0
        do i=1,nev_calc
            eigval(i)%val = eigval_all(ind(i))
            eigval(i)%prob = pev(ind(i))
            eigval(i)%sector = (ind(i)-1)/nev+1
            eigval(i)%level = mod(ind(i)-1,nev)+1

            Z = Z + pev(ind(i))
        enddo

        ! Normalizing the boltzman factor
        do i=1,nev_calc
            eigval(i)%prob = eigval(i)%prob / Z
        enddo

        if (node.eq.0) then
            write(6,*)
            write(6,*) "Obatined eigenvalues for all sectors."
            write(6,"(a,ES10.3,a,I5)") " Number of eigenvalues (with prob > ", &
                                       PROB_THRESHOLD,") = ", nev_calc
            write(6,*)
            write(6,"(a)") " Eigenvalue          Prob         Sector   Level"
            do i=1,nev_calc
                write(6,"(1x,ES16.5,4x,ES11.5,2I8)") eigval(i)%val,eigval(i)%prob,&
                                          eigval(i)%sector,eigval(i)%level
            enddo
            write(6,*)
        endif
    end subroutine ed_solve

    subroutine diag(basis, eigval, eigvec)
        ! include 'debug.h'
        include 'stat.h'
        type(basis_t), intent(in) :: basis
        double precision, intent(out) :: eigvec(basis%nloc,nev)
        double precision, intent(out) :: eigval(nev)

        integer  maxnloc,maxnev,maxncv,ldv
        parameter  (maxnloc=42000000,maxnev=40,maxncv=60,ldv=maxnloc)
        double precision :: v(ldv,maxncv), workl(maxncv*(maxncv+8)),        &
                    workd(3*maxnloc), d(maxncv,2), resid(maxnloc),  &
                    ax(maxnloc)
        logical :: select(maxncv)
        character, parameter :: bmat = 'I'
        character(len=2), parameter :: which = 'SA'
        integer :: iparam(11), ipntr(11), lworkl, info, ido, nconv, i, j, ncv, &
                   maxitr, mode, ishfts
        double precision :: sigma, tol
        double precision :: pdnorm2
        external :: pdnorm2, daxpy

        ! ndigit = -3
        ! logfil = 6
        ! msaupd = 3

        ncv = nev*2

        lworkl = ncv*(ncv+8)
        tol = 0.0

        ishfts = 1
        maxitr = 500
        mode   = 1

        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode

        info = 0
        ido = 0

        do
            call pdsaupd( comm, ido, bmat, basis%nloc, which, nev, tol, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if (ido .eq. -1 .or. ido .eq. 1) then
                call multiply_H(basis,workd(ipntr(1)),workd(ipntr(2)))
            else
                exit
            endif
        enddo

        if ( info .lt. 0 ) then
            if ( node.eq. 0 ) then
                print *, ' Error with pdsaupd, info = ', info
                print *, iparam(5)
                call die
            endif
        else
            call pdseupd(comm,.true.,'All', select,d,v,ldv,sigma, &
                bmat, basis%nloc, which, nev, tol, resid, ncv, v, ldv, &
                iparam, ipntr, workd, workl, lworkl, ierr )
            if ( ierr .ne. 0) then
                if ( node .eq. 0 ) then
                    print *, ' Error with pdseupd, info = ', ierr
                    call die
                endif
            else
                nconv =  iparam(5)
                ! do j=1, nconv
                !     call multiply_H(basis,v(:,j),ax)
                !     call daxpy(basis%nloc, -d(j,1), v(1,j), 1, ax, 1)
                !     d(j,2) = pdnorm2( comm, basis%nloc, ax, 1 )
                ! enddo
                ! call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                !     'Ritz values and direct residuals')
            end if
        endif

        eigval(:) = huge(1.0D0)
        do i=1,nconv
            eigvec(1:basis%nloc,i) = v(:,i)
            eigval(i) = d(i,1)
        enddo

    end subroutine diag

    subroutine boltzmann_factor(elv,nev,beta,p)

    implicit none

    integer:: i, k, nev
    double precision:: elv(nev), beta, p(nev), Z, en,betatmp

!      write(6,*) "This calculation assumes zero temperature!!!!"
!      betatmp = beta
!      beta = 10000*beta
    do i = 1, nev
       Z = 0.D0
       en = elv(i)
       do k = 1, nev
          if(-beta*(elv(k)-en).gt.10.D0) then
            p(i) = 0.D0
            goto 100
          else
            Z = Z + exp(-beta*(elv(k)-en))
          endif
       enddo
       p(i) = 1.D0/Z
100     continue
    enddo
!      beta = betatmp
    return
    end subroutine boltzmann_factor
end module ed_diag
