module ed_projection

    use matsubara_grid, only: omega, nwloc
    use ed_params, only: nbath
    use dmft_params, only: norb, mu, nw, nspin
    use mpi

    public :: project_to_impurity_model

    integer :: nx, iorb
    double complex, allocatable :: G0(:)

    private
contains

    ! obtain ek,vk by minimizing the difference between the cluster Weiss field and G0
    subroutine project_to_impurity_model(G0in,ek,vk)
        implicit none
        double complex, intent(in) :: G0in(nwloc,norb,nspin)
        double precision, intent(inout) :: ek(norb+nbath,nspin), vk(norb,nbath,nspin)

        ! local variables
        integer :: itr, i,j, ispin
        double precision, parameter :: tol = 1.0D-7 ! minimization tolerance 
        double precision :: x(1+2*nbath), diff

        ! For each (orb,spin),
        ! number of fitting parameters :
        !  - imp. lv : 1
        !  - bath lv : nbath
        !  - V_mk    : nbath
        ! Total : 1 + 2*nbath
        nx = 1 + 2*nbath

        allocate(G0(nwloc))

        if (master) then
            write(*,*)
            write(*,*) "Projecting to the impurity model..."
        endif

        do ispin=1,nspin
            do i=1,norb
                ! these following two variables are module-wide variables
                ! used in other subroutines
                iorb = i
                G0 = G0in(:,iorb,ispin)

                ! collect fitting parameters
                x(1) = ek(iorb,ispin)
                do j=1,nbath
                    x(1+j) = ek(norb+j,ispin)
                    x(1+nbath+j) = vk(iorb,j,ispin)
                enddo

                ! find x that minimizing the differences b/w G0 and G0cl
                call FRPRMN(x,nx,tol,itr,diff)

                ! recast x to ek,vk
                ek(iorb,ispin) = x(1)
                ek(norb+1:norb+nbath,ispin) = x(2:1+nbath)
                vk(iorb,:,ispin) = x(2+nbath:1+2*nbath)

                if (master) then
                    write(*,"(1x,a,I1,a,I1,a,I3,a,E12.4)") &
                        "(iorb,ispin)=(",iorb,",",ispin,") converged : itr=",&
                        itr,", diff=",diff
                endif

            enddo
        enddo


        if (master) then
            write(*,*) 
            write(*,*) "Projected impurity/bath Levels"
            do ispin=1,nspin
                write(*,*) "Spin ",ispin
                do i=1,norb+nbath
                    write(*,"(1x,A,I2,A,F12.6)") "ek(",i,") = ",ek(i,ispin)
                enddo
                write(*,*)
                write(*,*) "Impurity/Bath Hybridization"
                write(*,"(7x)", advance="no")
                do i=1,norb
                    write(*,"(4x,A3,I1,4x)",advance="no") "orb",i
                enddo
                write(*,*)
                do i=1,nbath
                    write(*,"(1x,a4,I2)",advance="no") "bath",i
                    do j=1,norb
                        write(*,"(F12.6)",advance="no") vk(j,i,ispin)
                    enddo
                    write(*,*)
                enddo
            enddo
            write(*,*)
        endif

        deallocate(G0)
    end subroutine project_to_impurity_model

    ! cluster hybridization function evaluation 
    double complex function delta_cl(iw,x)
        integer :: iw
        double precision :: x(nx)

        integer :: ibath, xidx

        delta_cl = cmplx(0.0d0,0.0d0)

        do ibath=1,nbath
            ! index for vk(iorb,ibath) in x
            xidx = 1+nbath+ibath
            delta_cl = delta_cl + x(xidx)*x(xidx)/(cmplx(0.0d0,omega(iw))-x(1+ibath))
        enddo

    end function delta_cl

    double precision function func(x)
        implicit none
        double precision :: x(nx)
        double precision :: func_loc, diff
        double complex :: gf

        integer :: iw, ibath, xidx

        func_loc = 0.0d0
        do iw=1,nwloc
            gf = cmplx(0.0d0,omega(iw))+mu-x(1)-delta_cl(iw,x)
            gf = 1/gf
            diff = abs(gf - G0(iw))

            ! weight 1/omega(iw)
            func_loc = func_loc + diff*diff/omega(iw)
        enddo

        call mpi_allreduce(func_loc,func,1,mpi_double_precision,mpi_sum,comm,mpierr)
    end function func

    ! overall factor for d/dxi |g0cl-g0|^2
    ! 2/(A^2 + B^2)^2
    double precision function fac0(A,B)
        double precision :: A,B
        fac0 = 2/(A**2+B**2)**2
    end function fac0

    ! coefficient of A'
    ! ( A^2- B^2 ) Re(G0(iw)) - A (1 + 2 B Im(G0(iw)))
    double precision function fac1(A,B,iw)
        double precision :: A,B
        integer :: iw
        fac1 = (A**2-B**2)*real(g0(iw))-A*(1+2*B*aimag(g0(iw)))
    end function fac1

    ! coefficient of B'
    ! B (-1 + 2 A Re(G0)) + ( A^2 - B^2 ) Im(G0)
    double precision function fac2(A,B,iw)
        double precision :: A,B
        integer :: iw
        fac2 = B*(-1+2*A*real(g0(iw))) - (A**2-B**2)*aimag(g0(iw))
    end function fac2

    subroutine dfunc(x,df)
        implicit none
        double precision :: x(nx),df(nx),df_loc(nx)

        double precision :: atmp, diff, A, B, P 
        double complex :: gf, dcl

        integer :: iw, ibath, ekidx, vkidx

        df_loc(1) = 0.0d0
        do iw=1,nwloc
            dcl = delta_cl(iw,x)
            A = mu - x(1) - real(dcl)
            B = omega(iw) - aimag(dcl)

            ! Ap = - delta(m,m')
            ! Bp = 0
            ! weight = 1/omega(iw)
            df_loc(1) = df_loc(1) - fac0(A,B)*fac1(A,B,iw)/omega(iw)
        enddo

        ! derivative wrt bath levels
        do ibath=1,nbath
            ekidx = 1+ibath
            vkidx = 1+nbath+ibath
            df_loc(ekidx) = 0.0d0
            do iw=1,nwloc
                ! Ap = 0
                ! Bp = -v(iorb,ibath)**2 * omega(iw) * 2 * ek(norb+ibath)/
                !       (omega(iw)**2+ek(norb+ibath)**2)**2
                P = -x(vkidx)**2 * omega(iw) * 2 * x(ekidx)/ &
                    (omega(iw)**2+x(ekidx)**2)**2

                ! weight 1/omega(iw)
                df_loc(ekidx) = df_loc(ekidx) + fac0(A,B)/omega(iw)*fac2(A,B,iw)*P
            enddo
        enddo

        ! derivatives wrt hybridization..
        do ibath=1,nbath
            vkidx = 1 + nbath + ibath
            ekidx = 1 + ibath
            df_loc(vkidx) = 0.0d0
            do iw=1,nwloc
                ! Ap = 0
                ! Bp = 2*x(vkidx)*omega(iw)/(omega(iw)**2+ek(norb+ibath)**2)
                P = 2*x(vkidx)*omega(iw)/(omega(iw)**2+x(ekidx)**2)
                df_loc(vkidx) = df_loc(vkidx)+fac0(A,B)/omega(iw)*fac2(A,B,iw)*P
            enddo
        enddo

        call mpi_allreduce(df_loc,df,nx,mpi_double_precision,mpi_sum,comm,mpierr)
    end subroutine dfunc

!=============================================================================================================================
!=============================================================================================================================
!     The following routines are from the numerical recipes.
!     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
!     is performed on a function func, using its gradient as calculated by a routine dfunc.
!     The convergence tolerance on the function value is input as ftol. Returned quantities are p 
!     (the location of the minimum), iter ( the # of iterations that were performed), and fret (the minimum value of function)
!     The routine linmin is called to perform line minimizations
!     parameter: nmax is the maximum anticipated value of n; itmax is the maximum allowd number of iterations;
!     eps is a small number of rectify special case of converging to exactly zero function value.
!=============================================================================================================================
!=============================================================================================================================
      SUBROUTINE FRPRMN(P,N,FTOL,ITER,FRET)
      PARAMETER (NMAX=50,ITMAX=200,EPS=1.E-10)
      DOUBLE PRECISION P(N),G(NMAX),H(NMAX),XI(NMAX)
      DOUBLE PRECISION FRET, FTOL
      FP=FUNC(P)
      CALL DFUNC(P,XI)
      DO 11 J=1,N
        G(J)=-XI(J)
        H(J)=G(J)
        XI(J)=H(J)
11    CONTINUE
      DO 14 ITS=1,ITMAX
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET)
        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN
        FP=FUNC(P)
        CALL DFUNC(P,XI)
        GG=0.
        DGG=0.
        DO 12 J=1,N
          GG=GG+G(J)**2
!         DGG=DGG+XI(J)**2
          DGG=DGG+(XI(J)+G(J))*XI(J)
12      CONTINUE
        IF(GG.EQ.0.)RETURN
        GAM=DGG/GG
        DO 13 J=1,N
          G(J)=-XI(J)
          H(J)=G(J)+GAM*H(J)
          XI(J)=H(J)
13      CONTINUE
14    CONTINUE
      STOP 'FRPR maximum iterations exceeded'
      RETURN
      END SUBROUTINE FRPRMN

      SUBROUTINE LINMIN(P,XI,N,FRET)
      PARAMETER (NMAX=50,TOL=1.E-4)
      ! EXTERNAL F1DIM
      DOUBLE PRECISION P(N),XI(N),FRET
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      NCOM=N
      DO 15 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
15    CONTINUE
      AX=0.
      XX=1.
      BX=2.
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      DO 16 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
16    CONTINUE
      RETURN
      END SUBROUTINE LINMIN

      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20)
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
17    IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 17
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 17
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 17
      ENDIF
      RETURN
      END SUBROUTINE MNBRAK

      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 18 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 23
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.P.GE.Q*(B-X)) GOTO 25
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 24
        ENDIF
25      IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
24      IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
18    CONTINUE
      STOP 'Brent exceed maximum iterations.'
23    XMIN=X
      BRENT=FX
      RETURN
      END FUNCTION BRENT

      FUNCTION F1DIM(X)
      PARAMETER (NMAX=50)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      DOUBLE PRECISION XT(NMAX)
      DO 30 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
30    CONTINUE
      F1DIM=FUNC(XT)
      RETURN
      END FUNCTION F1DIM
end module ed_projection

